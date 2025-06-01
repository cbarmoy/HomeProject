import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from collections import deque
import copy
import warnings

# Suppress deprecation warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

# Dash imports with error handling
try:
    import dash
    from dash import dcc, html
    from dash.dependencies import Input, Output, State
    import plotly.graph_objs as go
    from plotly.subplots import make_subplots
    import pandas as pd
    import dash_bootstrap_components as dbc
except ImportError as e:
    print(f"Error importing Dash dependencies: {e}")
    print("Please make sure you have installed all required packages:")
    print("pip install dash plotly pandas dash-bootstrap-components")
    raise

# ====================
# --- Stream Flow Field ---
# ====================

class StreamField:
    """Simplified flow field using streamlines and potential flow"""
    def __init__(self, width, height):
        self.width = width
        self.height = height
        
        # Adjusted flow parameters for better suspension
        self.base_flow = -0.5e-6  # Reduced downward flow
        self.drain_strength = -1e-6  # Reduced drain suction
        self.upward_flow_strength = 1e-6  # Added upward flow component
        self.drain_positions = [
            (0.3, 0.4),  # (start_x, end_x) for first drain
            (0.6, 0.7)   # (start_x, end_x) for second drain
        ]
        
    def get_velocity_at_point(self, x, y):
        """Calculate velocity at any point using superposition of flows"""
        # Initialize with gentler base flow
        vx, vy = 0.0, self.base_flow
        
        # Add upward flow component that increases with depth
        rel_y = y / self.height
        upward_flow = self.upward_flow_strength * (1 - rel_y)  # Stronger at bottom
        vy += upward_flow
        
        # Add circulation effects
        for start_x, end_x in self.drain_positions:
            # Convert to absolute coordinates
            drain_center_x = (start_x + end_x) / 2 * self.width
            drain_y = 0  # Drain at bottom
            
            # Calculate distance to drain center
            dx = x - drain_center_x
            dy = y - drain_y
            distance = np.hypot(dx, dy)
            
            if distance > 0:
                # Gentler drain effect
                strength = self.drain_strength / (1 + distance * 2e5)  # Reduced effect range
                vx += dx / distance * strength
                vy += dy / distance * strength
        
        return vx, vy

    def is_near_drain(self, x, y):
        """Check if a point is near any drain"""
        rel_x = x / self.width
        rel_y = y / self.height
        
        # Check if point is near bottom and within any drain x-range
        if rel_y < 0.08:  # Within 8% of bottom
            for start_x, end_x in self.drain_positions:
                if start_x <= rel_x <= end_x:
                    return True
        return False

# ====================
# --- Spatial Grid ---
# ====================

class SpatialGrid:
    """Spatial partitioning grid for efficient neighbor searches"""
    def __init__(self, width, height, cell_size):
        self.width = width
        self.height = height
        self.cell_size = cell_size
        
        # Calculate grid dimensions
        self.cols = int(np.ceil(width / cell_size))
        self.rows = int(np.ceil(height / cell_size))
        
        # Initialize grid
        self.grid = {}  # Dictionary for sparse storage
        
    def _get_cell_index(self, x, y):
        """Get grid cell indices for a position"""
        col = int(x / self.cell_size)
        row = int(y / self.cell_size)
        return min(max(row, 0), self.rows - 1), min(max(col, 0), self.cols - 1)
    
    def update_particle(self, particle, particle_id):
        """Update particle position in grid"""
        cell_index = self._get_cell_index(particle.x, particle.y)
        if cell_index not in self.grid:
            self.grid[cell_index] = set()
        self.grid[cell_index].add(particle_id)
    
    def clear(self):
        """Clear all particles from grid"""
        self.grid.clear()
    
    def get_nearby_indices(self, x, y, radius):
        """Get indices of particles within radius of position"""
        # Calculate cell range to check
        start_row, start_col = self._get_cell_index(x - radius, y - radius)
        end_row, end_col = self._get_cell_index(x + radius, y + radius)
        
        nearby = set()
        for row in range(start_row, end_row + 1):
            for col in range(start_col, end_col + 1):
                if (row, col) in self.grid:
                    nearby.update(self.grid[(row, col)])
        return nearby

def create_spatial_grids(space_size):
    """Create spatial grids for different particle types"""
    # Use different cell sizes based on typical search radii
    bacteria_grid = SpatialGrid(space_size, space_size, cell_size=5e-6)  # Larger cells for bacteria
    glucose_grid = SpatialGrid(space_size, space_size, cell_size=2e-6)   # Smaller cells for glucose
    metabolite_grid = SpatialGrid(space_size, space_size, cell_size=3e-6)  # Medium cells for metabolites
    return bacteria_grid, glucose_grid, metabolite_grid

def update_spatial_grids(bacteria_grid, glucose_grid, metabolite_grid, 
                        ecoli_list, glucose_list, metabolite_list):
    """Update all spatial grids with current particle positions"""
    bacteria_grid.clear()
    glucose_grid.clear()
    metabolite_grid.clear()
    
    for i, bacteria in enumerate(ecoli_list):
        bacteria_grid.update_particle(bacteria, i)
    
    for i, glucose in enumerate(glucose_list):
        glucose_grid.update_particle(glucose, i)
    
    for i, metabolite in enumerate(metabolite_list):
        metabolite_grid.update_particle(metabolite, i)

# ====================
# --- Classes ---
# ====================

@dataclass
class DrainSystem:
    """System for collecting and recycling particles through filtration"""
    space_size: float
    collected_metabolites: int = 0
    drain_positions: list = field(default_factory=lambda: [(0.3, 0.4), (0.6, 0.7)])
    drain_height: float = 0.08
    collection_radius: float = 0.02
    recycling_queue: list = field(default_factory=list)  # Queue for recycling
    recycling_delay: int = 25  # Steps to wait before recycling (5 seconds at 200ms interval)
    
    def check_collection(self, particle, particle_type):
        """Check if a particle should be collected by the drain"""
        rel_x = particle.x / self.space_size
        rel_y = particle.y / self.space_size
        
        # Check if within drain height
        if rel_y < self.drain_height:
            for start_x, end_x in self.drain_positions:
                drain_center_x = (start_x + end_x) / 2
                distance_to_center = abs(rel_x - drain_center_x)
                
                if distance_to_center <= (end_x - start_x) / 2 + self.collection_radius:
                    collection_probability = 1 - (rel_y / self.drain_height)
                    if np.random.random() < collection_probability:
                        if particle_type == "metabolite":
                            self.collected_metabolites += 1
                        else:
                            # Add to recycling queue with timestamp and type
                            self.recycling_queue.append({
                                'particle': copy.deepcopy(particle),
                                'type': particle_type,
                                'collected_at': len(self.recycling_queue)
                            })
                        return True
        return False
    
    def process_recycling(self, current_step):
        """Process particles ready for recycling"""
        if not self.recycling_queue:
            return [], []  # Return only two empty lists when queue is empty
        
        new_bacteria = []
        new_glucose = []
        ready_indices = []
        
        for i, item in enumerate(self.recycling_queue):
            if current_step - item['collected_at'] >= self.recycling_delay:
                particle = item['particle']
                
                # Reset position to top of reactor with random x position
                particle.x = np.random.uniform(0.2 * self.space_size, 0.8 * self.space_size)
                particle.y = np.random.uniform(0.85 * self.space_size, 0.95 * self.space_size)
                
                # Reset momentum
                if hasattr(particle, 'vertical_momentum'):
                    particle.vertical_momentum = 0
                    particle.horizontal_momentum = np.random.normal(0, 0.5e-6)
                
                if item['type'] == 'bacteria':
                    new_bacteria.append(particle)
                elif item['type'] == 'glucose':
                    new_glucose.append(particle)
                
                ready_indices.append(i)
        
        # Remove processed items from queue
        for i in sorted(ready_indices, reverse=True):
            self.recycling_queue.pop(i)
        
        return new_bacteria, new_glucose

@dataclass
class Metabolite:
    x: float
    y: float
    type: str = "generic"
    diameter: float = 0.0005e-6
    diffusion_rate: float = 2.5e-7
    decay_rate: float = 0.003
    concentration: float = 0.3
    toxicity_threshold: float = 8.0
    max_toxicity: float = 0.3
    gravity: float = 6e-6  # Added gravity effect for metabolites
    
    def calculate_toxicity(self, local_concentration):
        """Calculate toxic effect based on local concentration"""
        if local_concentration > self.toxicity_threshold:
            # Even more gradual toxicity increase
            relative_conc = (local_concentration - self.toxicity_threshold) / (self.toxicity_threshold * 3)
            toxicity = self.max_toxicity * (1 / (1 + np.exp(-relative_conc + 4)))
            return toxicity * 0.4
        return 0.0
    
    def diffuse(self, diffusion_strength, space_size, fluid_field=None):
        """Diffuse metabolite through the medium with gravity"""
        # Random diffusion with higher rate
        self.x += np.random.normal(0, self.diffusion_rate * diffusion_strength)
        self.y += np.random.normal(0, self.diffusion_rate * diffusion_strength)
        
        # Apply gravity
        self.y -= self.gravity * diffusion_strength * 2
        
        # Apply fluid velocity if fluid field exists
        if fluid_field is not None:
            u, v = fluid_field.get_velocity_at_point(self.x, self.y)
            self.x += u * diffusion_strength * 150  # Increased fluid effect
            self.y += v * diffusion_strength * 150
        
        # Ensure staying within bounds
        self.x = np.clip(self.x, 0, space_size)
        self.y = np.clip(self.y, 0, space_size)
        
        # Increase concentration decay over time
        self.concentration *= (1 - self.decay_rate * 1.2)

    def __deepcopy__(self, memo):
        return Metabolite(
            x=self.x,
            y=self.y,
            type=self.type,
            diameter=self.diameter,
            diffusion_rate=self.diffusion_rate,
            decay_rate=self.decay_rate,
            concentration=self.concentration,
            toxicity_threshold=self.toxicity_threshold,
            max_toxicity=self.max_toxicity,
            gravity=self.gravity
        )

@dataclass
class Glucose:
    x: float
    y: float
    diameter: float = 0.0009e-6
    gravity: float = 0.005e-6  # Increased gravity for glucose
    vertical_momentum: float = 0.0
    horizontal_momentum: float = 0.0

    def diffuse(self, diffusion_strength, space_size, fluid_field=None):
        # Keep existing diffusion parameters but increase gravity effect
        base_diffusion = 2.5e-6
        momentum_decay = 0.95
        
        # Update momentum with random forces
        self.horizontal_momentum = (self.horizontal_momentum * momentum_decay + 
                                  np.random.normal(0, base_diffusion))
        self.vertical_momentum = (self.vertical_momentum * momentum_decay + 
                                np.random.normal(-0.1e-6, base_diffusion))  # Slight downward bias
        
        # Apply momentum to position
        self.x += self.horizontal_momentum * diffusion_strength
        self.y += self.vertical_momentum * diffusion_strength
        
        # Increased gravity effect
        gravity_effect = self.gravity * 0.8  # Much stronger gravity
        self.y -= gravity_effect * diffusion_strength
        
        # Apply fluid effects
        if fluid_field is not None:
            u, v = fluid_field.get_velocity_at_point(self.x, self.y)
            fluid_factor = 8.0
            
            self.horizontal_momentum += u * diffusion_strength * fluid_factor
            self.vertical_momentum += v * diffusion_strength * fluid_factor * 0.3
        
        # Boundary interactions
        bounce_factor = 0.7
        
        if self.x <= 0:
            self.x = 0
            self.horizontal_momentum = abs(self.horizontal_momentum) * bounce_factor
        elif self.x >= space_size:
            self.x = space_size
            self.horizontal_momentum = -abs(self.horizontal_momentum) * bounce_factor
            
        if self.y <= 0:
            self.y = 0
            self.vertical_momentum = abs(self.vertical_momentum) * bounce_factor
        elif self.y >= space_size:
            self.y = space_size
            self.vertical_momentum = -abs(self.vertical_momentum) * bounce_factor

    def __deepcopy__(self, memo):
        new_glucose = Glucose(x=self.x, y=self.y, diameter=self.diameter, gravity=self.gravity)
        new_glucose.vertical_momentum = self.vertical_momentum
        new_glucose.horizontal_momentum = self.horizontal_momentum
        return new_glucose

@dataclass
class Ecoli:
    x: float
    y: float
    ATP: float = 200
    width: float = 2e-6
    height: float = 1e-6
    produced_metabolite: int = 0
    gravity: float = 0.5e-6  # Significantly reduced gravity effect
    glucose_sensitivity: float = 0.5
    base_growth_rate: float = 0.12
    metabolite_excretion_threshold: float = 120
    metabolite_excretion_cost: float = 4
    metabolite_sensing_radius: float = 5e-6
    
    def calculate_growth_probability(self, glucose_list, glucose_grid, sensing_radius=5e-6):
        """Calculate growth probability using spatial grid for efficiency"""
        nearby_indices = glucose_grid.get_nearby_indices(self.x, self.y, sensing_radius)
        # Filter out invalid indices
        valid_indices = [i for i in nearby_indices if i < len(glucose_list)]
        local_glucose = [glucose_list[i] for i in valid_indices]
        
        # Calculate probability based on local glucose count
        glucose_count = len(local_glucose)
        probability = self.base_growth_rate * (glucose_count / (1/self.glucose_sensitivity + glucose_count))
        
        return probability, local_glucose
    
    def calculate_local_metabolite_concentration(self, metabolite_list, metabolite_grid):
        """Calculate local metabolite concentration using spatial grid"""
        nearby_indices = metabolite_grid.get_nearby_indices(
            self.x, self.y, self.metabolite_sensing_radius
        )
        
        local_concentration = 0
        for idx in nearby_indices:
            # Skip invalid indices
            if idx >= len(metabolite_list):
                continue
                
            metabolite = metabolite_list[idx]
            distance = np.hypot(self.x - metabolite.x, self.y - metabolite.y)
            if distance <= self.metabolite_sensing_radius:
                weight = (1 - (distance / self.metabolite_sensing_radius)) ** 3
                local_concentration += metabolite.concentration * weight * 0.5
        
        return local_concentration
    
    def apply_toxic_effects(self, metabolite_list, metabolite_grid):
        """Apply toxic effects from high metabolite concentration"""
        local_concentration = self.calculate_local_metabolite_concentration(metabolite_list, metabolite_grid)
        
        # Get maximum toxicity effect from all metabolites
        max_toxic_effect = 0
        for metabolite in metabolite_list:
            toxic_effect = metabolite.calculate_toxicity(local_concentration)
            max_toxic_effect = max(max_toxic_effect, toxic_effect)
        
        if max_toxic_effect > 0:
            # Even more gradual ATP reduction
            atp_reduction = self.ATP * max_toxic_effect * 0.2
            self.ATP -= atp_reduction
            return atp_reduction
        return 0

    def excrete_metabolite(self, space_size):
        """Excrete metabolite if enough ATP is available"""
        if self.ATP >= self.metabolite_excretion_threshold:
            # Create new metabolite near the bacteria
            offset_x = np.random.uniform(-self.width/2, self.width/2)
            offset_y = np.random.uniform(-self.height/2, self.height/2)
            
            metabolite = Metabolite(
                x=np.clip(self.x + offset_x, 0, space_size),
                y=np.clip(self.y + offset_y, 0, space_size),
                type="primary_metabolite"
            )
            
            # Deduct ATP cost
            self.ATP -= self.metabolite_excretion_cost
            self.produced_metabolite += 1
            
            return metabolite
        return None

    def move(self, glucose_list, glucose_grid, step_size, noise, space_size, detection_radius, fluid_field=None):
        """Stochastic movement with chemotaxis"""
        # Get nearby glucose for chemotaxis
        _, nearby_glucose = self.calculate_growth_probability(glucose_list, glucose_grid, detection_radius)
        
        # Add stochastic component (noise) - Moved up here
        ndx, ndy = np.random.uniform(-1, 1), np.random.uniform(-1, 1)
        
        if nearby_glucose:
            # Move towards nearest glucose (chemotaxis)
            nearest = min(nearby_glucose, key=lambda g: (self.x - g.x) ** 2 + (self.y - g.y) ** 2)
            dx, dy = nearest.x - self.x, nearest.y - self.y
            norm = np.hypot(dx, dy)
            if norm > 0:
                dx /= norm
                dy /= norm
        else:
            # Random movement if no glucose detected
            dx, dy = np.random.uniform(-1, 1), np.random.uniform(-1, 1)
            norm = np.hypot(dx, dy)
            if norm > 0:
                dx /= norm
                dy /= norm
        
        # Combine directed movement, noise, and very gentle gravity
        final_dx = step_size * ((1 - noise) * dx + noise * ndx)
        final_dy = step_size * ((1 - noise) * dy + noise * ndy) - self.gravity * step_size * 2  # Reduced gravity effect
        
        # Apply fluid velocity if fluid field exists
        if fluid_field is not None:
            u, v = fluid_field.get_velocity_at_point(self.x, self.y)
            final_dx += u * step_size * 15
            final_dy += v * step_size * 10  # Reduced vertical fluid effect
        
        # Update position with bounds checking
        self.x = np.clip(self.x + final_dx, 0, space_size)
        self.y = np.clip(self.y + final_dy, 0, space_size)

        # ATP costs
        # Reduced cost for moving against gravity
        if final_dy > 0:  # Moving upward
            self.ATP -= abs(final_dy) * 0.5  # Reduced upward movement cost
        
        # Base movement cost
        movement_cost = np.hypot(final_dx, final_dy) * 0.5
        self.ATP -= movement_cost

    def consume(self, glucose, atp_gain):
        """Discrete glucose consumption event"""
        threshold = (glucose.diameter / 2) + max(self.width, self.height) / 2
        if np.hypot(self.x - glucose.x, self.y - glucose.y) <= threshold:
            self.ATP += atp_gain * 0.8  # Reduced ATP gain from glucose
            self.produced_metabolite += 1
            return True
        return False

    def can_divide(self, threshold):
        """Check if cell has enough ATP to divide"""
        return self.ATP >= threshold

    def divide(self, space_size):
        """Discrete cell division event"""
        self.ATP /= 2
        # New bacteria appears slightly to the side and above/below
        offset_x = np.random.uniform(-1e-6, 1e-6)
        offset_y = np.random.uniform(-1e-6, 1e-6)
        return Ecoli(
            x=np.clip(self.x + offset_x, 0, space_size),
            y=np.clip(self.y + offset_y, 0, space_size),
            ATP=self.ATP,
            gravity=self.gravity,
            glucose_sensitivity=self.glucose_sensitivity,
            base_growth_rate=self.base_growth_rate,
            metabolite_excretion_threshold=self.metabolite_excretion_threshold,
            metabolite_excretion_cost=self.metabolite_excretion_cost
        )

    def is_alive(self, threshold):
        """Check if cell has enough ATP to stay alive"""
        return self.ATP >= threshold

    def __deepcopy__(self, memo):
        return Ecoli(
            x=self.x,
            y=self.y,
            ATP=self.ATP,
            width=self.width,
            height=self.height,
            produced_metabolite=self.produced_metabolite,
            gravity=self.gravity,
            glucose_sensitivity=self.glucose_sensitivity,
            base_growth_rate=self.base_growth_rate,
            metabolite_excretion_threshold=self.metabolite_excretion_threshold,
            metabolite_excretion_cost=self.metabolite_excretion_cost
        )

# ====================
# --- Gillespie algorithm ---
# ====================

def gillespie_step(ecoli, glucose_list, glucose_grid, space_size):
    """
    Gillespie algorithm step for discrete event simulation
    Events: division, death, or idle (no change)
    """
    # Calculate growth probability based on local glucose
    growth_probability, _ = ecoli.calculate_growth_probability(glucose_list, glucose_grid)
    
    # ATP decay over time (bacteria consume energy to stay alive)
    ecoli.ATP -= 0.5  # Reduced base metabolism cost
    
    # Calculate death probability based on ATP level
    death_probability = 0.005  # Reduced chance of random death
    if ecoli.ATP < 50:  # Lower ATP threshold for death
        death_probability = 0.1  # Reduced base death rate when ATP is low
        if ecoli.ATP < 25:
            death_probability = 0.2  # Reduced death rate when ATP is very low
        if ecoli.ATP < 10:
            death_probability = 0.5  # Reduced death rate when ATP is critically low
    
    # Define possible events and their rates
    rates = {
        'divide': growth_probability if ecoli.can_divide(200) else 0,  # Increased division probability
        'die': death_probability,
        'idle': 0.9  # Increased idle rate for more stability
    }
    
    # Calculate total rate
    total_rate = sum(rates.values())
    if total_rate == 0:
        return 'idle', None
    
    # Generate random numbers for event selection
    rand = np.random.uniform(0, total_rate)
    
    # Select event based on probability
    cumulative = 0
    for action, rate in rates.items():
        cumulative += rate
        if rand < cumulative:
            if action == 'divide':
                return action, ecoli.divide(space_size)
            return action, None
            
    return 'idle', None

# ====================
# --- Simulation ---
# ====================

def add_glucose_from_top(glucose_list, n_glucose, space_size):
    """
    Add new glucose molecules with improved floating behavior
    """
    new_glucose = []
    
    # Define zones with more emphasis on upper regions
    zones = [
        (0.90, 1.0, 0.3),     # Top zone: 30% of new glucose
        (0.75, 0.90, 0.3),    # Upper zone: 30% of new glucose
        (0.60, 0.75, 0.25),   # Middle zone: 25% of new glucose
        (0.45, 0.60, 0.15),   # Lower zone: 15% of new glucose
    ]
    
    for zone_top, zone_bottom, proportion in zones:
        n_zone_glucose = int(n_glucose * proportion)
        for _ in range(n_zone_glucose):
            # Wide horizontal distribution
            x_position = np.random.uniform(0.1 * space_size, 0.9 * space_size)
            
            # Vertical position within zone
            y_position = np.random.uniform(
                zone_bottom * space_size,
                zone_top * space_size
            )
            
            # Create glucose with initial upward momentum
            glucose = Glucose(x_position, y_position)
            glucose.vertical_momentum = np.random.uniform(0.5e-6, 1.5e-6)  # Initial upward momentum
            glucose.horizontal_momentum = np.random.normal(0, 1e-6)  # Random horizontal momentum
            
            new_glucose.append(glucose)
    
    glucose_list.extend(new_glucose)
    return glucose_list

def create_ecoli_and_glucose(n_ecoli, n_glucose, space_size):
    ecoli_list = [Ecoli(np.random.uniform(0, space_size), np.random.uniform(0, space_size))
                  for _ in range(n_ecoli)]
    glucose_list = [Glucose(np.random.uniform(0, space_size), np.random.uniform(0, space_size))
                    for _ in range(n_glucose)]
    return ecoli_list, glucose_list

def simulate(n_ecoli=10, n_glucose=300, space_size=20e-6, max_steps=250, 
            glucose_feed_interval=20, glucose_feed_amount=50):
    """
    Simulate bacterial growth with simplified fluid dynamics and discrete events
    """
    # Initialize spatial grids
    bacteria_grid, glucose_grid, metabolite_grid = create_spatial_grids(space_size)
    
    # Initialize stream field for fluid flow
    fluid = StreamField(width=space_size, height=space_size)
    
    # Initialize drain system
    drain_system = DrainSystem(space_size=space_size)
    
    ecoli_list, glucose_list = create_ecoli_and_glucose(n_ecoli, n_glucose, space_size)
    metabolite_list = []
    states = []

    # Initialize statistics tracking
    total_toxic_deaths = 0
    total_atp_lost_to_toxicity = 0

    # Print initial setup
    print("\nStarting simulation...")
    print(f"Initial bacteria: {n_ecoli}")
    print(f"Initial glucose: {n_glucose}")
    print(f"Total steps: {max_steps}\n")

    # Dynamic glucose feed adjustment based on bacterial population
    base_feed_amount = glucose_feed_amount
    target_glucose = n_glucose

    for step in range(max_steps):
        # Calculate and display progress
        progress = (step / max_steps) * 100
        print(f"\rProgress: {progress:3.1f}% | Step: {step}/{max_steps} | "
              f"Bacteria: {len(ecoli_list)} | Glucose: {len(glucose_list)} | "
              f"Metabolites: {len(metabolite_list)} | Collected: {drain_system.collected_metabolites}", 
              end="", flush=True)
        
        # Process recycling
        new_bacteria, new_glucose = drain_system.process_recycling(step)
        ecoli_list.extend(new_bacteria)
        glucose_list.extend(new_glucose)
        
        # Update spatial grids
        update_spatial_grids(bacteria_grid, glucose_grid, metabolite_grid,
                           ecoli_list, glucose_list, metabolite_list)
        
        # Process metabolites with drain collection
        if metabolite_list:
            new_metabolites = []
            for metabolite in metabolite_list:
                if not drain_system.check_collection(metabolite, "metabolite"):
                    # Apply diffusion before adding to new list
                    metabolite.diffuse(diffusion_strength=1.0, space_size=space_size, fluid_field=fluid)
                    new_metabolites.append(metabolite)
            metabolite_list = new_metabolites
        
        # Process glucose with drain collection
        if glucose_list:
            new_glucose = []
            for glucose in glucose_list:
                if not drain_system.check_collection(glucose, "glucose"):
                    # Apply diffusion before adding to new list
                    glucose.diffuse(diffusion_strength=1.0, space_size=space_size, fluid_field=fluid)
                    new_glucose.append(glucose)
            glucose_list = new_glucose
        
        # Process bacteria with drain collection
        new_bact = []
        survivors = []
        for bacteria in ecoli_list:
            if not drain_system.check_collection(bacteria, "bacteria"):
                # Move bacteria
                bacteria.move(glucose_list, glucose_grid, step_size=1e-6, noise=0.5,
                            space_size=space_size, detection_radius=1e-5,
                            fluid_field=fluid)
                
                # Handle metabolite excretion
                new_metabolite = bacteria.excrete_metabolite(space_size)
                if new_metabolite:
                    metabolite_list.append(new_metabolite)
                
                # Handle division and death
                action, result = gillespie_step(bacteria, glucose_list, glucose_grid, space_size)
                if action == 'divide' and result:
                    new_bact.append(result)
                    survivors.append(bacteria)
                elif action != 'die':
                    survivors.append(bacteria)

        ecoli_list = survivors + new_bact
        
        # Dynamic glucose feeding with improved distribution
        if step > 0 and step % glucose_feed_interval == 0:
            population_factor = min(len(ecoli_list) / max(n_ecoli, 1), 2.0)
            current_glucose = len(glucose_list)
            
            if current_glucose < target_glucose * 0.1:
                glucose_factor = 2.0
            else:
                glucose_factor = max(0.5, min(1.5, target_glucose / max(current_glucose, 1)))
            
            adjusted_feed_amount = max(
                int(base_feed_amount * population_factor * glucose_factor),
                base_feed_amount // 2
            )
            glucose_list = add_glucose_from_top(glucose_list, adjusted_feed_amount, space_size)
        
        # Save state for visualization
        states.append((copy.deepcopy(ecoli_list), 
                      copy.deepcopy(glucose_list), 
                      copy.deepcopy(metabolite_list)))
        
        # Check for simulation end condition
        if not ecoli_list:
            print(f"\n\nSimulation stopped at step {step + 1}: no bacteria remaining.")
            break

    return states, space_size, drain_system

# ====================
# --- Data Export ---
# ====================

def export_detailed_data(states, filename="simulation_data.csv"):
    """
    Export detailed simulation data to CSV with comprehensive statistics
    
    The CSV includes the following categories of metrics:
    - Temporal: Step number and timestamp
    - Population: Counts of bacteria, glucose, and metabolites
    - ATP Metrics: Average, min, max, total ATP and standard deviation
    - Metabolite Metrics: Production, concentration, and collection statistics
    - Spatial Metrics: Average positions and distributions
    - Efficiency Metrics: Resource utilization and growth rates
    """
    import csv
    from datetime import datetime
    
    # Generate timestamp for the filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"simulation_data_{timestamp}.csv"
    
    statistics = []
    for step, (ecolis, glucoses, metabolites) in enumerate(states):
        # Skip if no bacteria (prevents division by zero)
        if not ecolis:
            continue
            
        # Basic counts
        n_bacteria = len(ecolis)
        n_glucose = len(glucoses)
        n_metabolites = len(metabolites)
        
        # ATP statistics
        atp_values = [e.ATP for e in ecolis]
        avg_atp = np.mean(atp_values)
        min_atp = min(atp_values)
        max_atp = max(atp_values)
        total_atp = sum(atp_values)
        atp_std = np.std(atp_values)
        
        # Metabolite statistics
        metabolite_production = [e.produced_metabolite for e in ecolis]
        total_metabolites_produced = sum(metabolite_production)
        avg_metabolites_per_bacteria = np.mean(metabolite_production)
        metabolite_concentration = [m.concentration for m in metabolites] if metabolites else [0]
        avg_metabolite_concentration = np.mean(metabolite_concentration)
        
        # Position statistics
        x_positions = [e.x for e in ecolis]
        y_positions = [e.y for e in ecolis]
        avg_bacteria_x = np.mean(x_positions)
        avg_bacteria_y = np.mean(y_positions)
        bacteria_spatial_std = np.std(x_positions) + np.std(y_positions)
        
        glucose_x = [g.x for g in glucoses]
        glucose_y = [g.y for g in glucoses]
        avg_glucose_x = np.mean(glucose_x) if glucose_x else 0
        avg_glucose_y = np.mean(glucose_y) if glucose_y else 0
        glucose_spatial_std = np.std(glucose_x) + np.std(glucose_y) if glucose_x else 0
        
        # Population dynamics
        prev_bacteria_count = len(states[step-1][0]) if step > 0 else n_bacteria
        growth_rate = (n_bacteria / prev_bacteria_count) if prev_bacteria_count > 0 else 1.0
        
        # Resource efficiency
        glucose_per_bacteria = n_glucose / n_bacteria if n_bacteria > 0 else 0
        atp_per_glucose = total_atp / n_glucose if n_glucose > 0 else 0
        metabolites_per_glucose = n_metabolites / n_glucose if n_glucose > 0 else 0
        
        stats = {
            # Temporal metrics
            'Step': step,
            'Time_stamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f"),
            
            # Population metrics
            'Bacteria_Count': n_bacteria,
            'Glucose_Count': n_glucose,
            'Metabolite_Count': n_metabolites,
            'Growth_Rate': f"{growth_rate:.3f}",
            
            # ATP metrics
            'Average_ATP': f"{avg_atp:.2f}",
            'Min_ATP': f"{min_atp:.2f}",
            'Max_ATP': f"{max_atp:.2f}",
            'Total_ATP': f"{total_atp:.2f}",
            'ATP_StdDev': f"{atp_std:.2f}",
            
            # Metabolite metrics
            'Total_Metabolites_Produced': total_metabolites_produced,
            'Avg_Metabolites_Per_Bacteria': f"{avg_metabolites_per_bacteria:.2f}",
            'Avg_Metabolite_Concentration': f"{avg_metabolite_concentration:.3f}",
            
            # Spatial distribution
            'Avg_Bacteria_X': f"{avg_bacteria_x:.6f}",
            'Avg_Bacteria_Y': f"{avg_bacteria_y:.6f}",
            'Bacteria_Spatial_Std': f"{bacteria_spatial_std:.6f}",
            'Avg_Glucose_X': f"{avg_glucose_x:.6f}",
            'Avg_Glucose_Y': f"{avg_glucose_y:.6f}",
            'Glucose_Spatial_Std': f"{glucose_spatial_std:.6f}",
            
            # Efficiency metrics
            'Glucose_Per_Bacteria': f"{glucose_per_bacteria:.2f}",
            'ATP_Per_Glucose': f"{atp_per_glucose:.2f}",
            'Metabolites_Per_Glucose': f"{metabolites_per_glucose:.2f}",
            'ATP_Per_Bacteria': f"{total_atp/n_bacteria:.2f}",
            'Metabolites_Per_Bacteria': f"{total_metabolites_produced/n_bacteria:.2f}"
        }
        statistics.append(stats)
    
    # Write to CSV with metrics grouped by category
    with open(filename, "w", newline="", encoding='utf-8') as csvfile:
        fieldnames = [
            # Temporal metrics
            'Step', 'Time_stamp',
            # Population metrics
            'Bacteria_Count', 'Glucose_Count', 'Metabolite_Count', 'Growth_Rate',
            # ATP metrics
            'Average_ATP', 'Min_ATP', 'Max_ATP', 'Total_ATP', 'ATP_StdDev',
            # Metabolite metrics
            'Total_Metabolites_Produced', 'Avg_Metabolites_Per_Bacteria', 'Avg_Metabolite_Concentration',
            # Spatial distribution
            'Avg_Bacteria_X', 'Avg_Bacteria_Y', 'Bacteria_Spatial_Std',
            'Avg_Glucose_X', 'Avg_Glucose_Y', 'Glucose_Spatial_Std',
            # Efficiency metrics
            'Glucose_Per_Bacteria', 'ATP_Per_Glucose', 'Metabolites_Per_Glucose',
            'ATP_Per_Bacteria', 'Metabolites_Per_Bacteria'
        ]
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for stats in statistics:
            writer.writerow(stats)
    
    print(f"\nDetailed simulation data exported to: {filename}")
    return filename


def export_summary_image(states, filename="simulation_summary.png"):
    """
    Generate and export a high-quality summary image of the simulation results.

    Parameters:
        states (list): List of simulation states over time.
        filename (str): Output image file name (PNG).
    """
    if not states:
        print("No simulation data to export.")
        return

    steps = list(range(len(states)))
    bacteria_counts = [len(s[0]) for s in states]
    glucose_counts = [len(s[1]) for s in states]
    metabolite_counts = [len(s[2]) for s in states]
    avg_atp = [np.mean([e.ATP for e in s[0]]) if s[0] else 0 for s in states]
    total_metabolites = [sum(e.produced_metabolite for e in s[0]) for s in states]

    plt.figure(figsize=(14, 10))

    # Subplot 1: Population dynamics
    plt.subplot(2, 2, 1)
    plt.plot(steps, bacteria_counts, label="E. coli", color='green')
    plt.plot(steps, glucose_counts, label="Glucose", color='blue')
    plt.plot(steps, metabolite_counts, label="Metabolites", color='purple')
    plt.title("Population Dynamics Over Time", fontsize=14)
    plt.xlabel("Simulation Step")
    plt.ylabel("Count")
    plt.legend()
    plt.grid(True)

    # Subplot 2: Average ATP
    plt.subplot(2, 2, 2)
    plt.plot(steps, avg_atp, label="Avg ATP", color='orange')
    plt.title("Average ATP Level in E. coli", fontsize=14)
    plt.xlabel("Simulation Step")
    plt.ylabel("ATP Units")
    plt.grid(True)

    # Subplot 3: Total metabolites produced
    plt.subplot(2, 2, 3)
    plt.plot(steps, total_metabolites, label="Cumulative Metabolites", color='teal')
    plt.title("Total Metabolites Produced", fontsize=14)
    plt.xlabel("Simulation Step")
    plt.ylabel("Metabolites")
    plt.grid(True)

    # Subplot 4: Ratios
    plt.subplot(2, 2, 4)
    ratio = [g / b if b > 0 else 0 for g, b in zip(glucose_counts, bacteria_counts)]
    plt.plot(steps, ratio, label="Glucose per Bacteria", color='red')
    plt.title("Glucose per Bacteria Ratio", fontsize=14)
    plt.xlabel("Simulation Step")
    plt.ylabel("Ratio")
    plt.grid(True)

    plt.suptitle("Simulation Summary", fontsize=18)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Export as high-resolution PNG
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"Summary image exported as {filename}")


# ====================
# --- Dash Visualization ---
# ====================

class SimulationDashboard:
    def __init__(self, states, space_size, fluid_field, drain_system):
        self.states = states
        self.space_size = space_size
        self.fluid_field = fluid_field
        self.drain_system = drain_system
        self.current_step = 0
        self.max_steps = len(states)
        self.reset_clicks = 0
        self.play_clicks = 0
        
        # Initialize data storage
        self.history_length = 100
        self.reset_histories()
        
        # Create Dash app with DARKLY theme
        external_stylesheets = [dbc.themes.DARKLY]
        self.app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
        
        # Store initial state
        self.initial_state = {
            'current_step': 0,
            'histories': {
                'time': list(self.time_history),
                'metabolites': list(self.metabolites_history),
                'population': list(self.population_history),
                'atp': list(self.atp_history),
                'excreted_metabolites': list(self.excreted_metabolites_history),
                'toxic_concentration': list(self.toxic_concentration_history)
            }
        }
        
        self.setup_layout()
        self.setup_callbacks()

    def reset_histories(self):
        """Reset all history deques"""
        self.metabolites_history = deque(maxlen=self.history_length)
        self.population_history = deque(maxlen=self.history_length)
        self.atp_history = deque(maxlen=self.history_length)
        self.time_history = deque(maxlen=self.history_length)
        self.excreted_metabolites_history = deque(maxlen=self.history_length)
        self.toxic_concentration_history = deque(maxlen=self.history_length)

    def calculate_average_toxicity(self, ecoli_list, metabolite_list):
        """Calculate average toxic concentration experienced by bacteria"""
        if not ecoli_list:
            return 0
        
        # Create a temporary metabolite grid for the calculation
        metabolite_grid = SpatialGrid(self.space_size, self.space_size, cell_size=3e-6)
        for i, metabolite in enumerate(metabolite_list):
            metabolite_grid.update_particle(metabolite, i)
        
        total_concentration = sum(
            e.calculate_local_metabolite_concentration(metabolite_list, metabolite_grid)
            for e in ecoli_list
        )
        return total_concentration / len(ecoli_list)

    def run(self, debug=False, port=8050):
        """Run the dashboard server"""
        if hasattr(self, 'app'):
            print("\n" + "="*50)
            print("Dashboard is running at: http://127.0.0.1:8050")
            print("Open this link in your web browser to view the simulation")
            print("="*50 + "\n")
            
            # Suppress Flask logging
            import logging
            log = logging.getLogger('werkzeug')
            log.setLevel(logging.ERROR)
            
            # Run the server with minimal output
            self.app.run(debug=debug, port=port, host='127.0.0.1', use_reloader=False)
        else:
            print("Error: Dashboard app not properly initialized")

    def setup_layout(self):
        """Setup the dashboard layout"""
        self.app.layout = dbc.Container([
            dbc.Row([
                dbc.Col([
                    html.H1("Bacterial Growth Simulation",
                           className="text-center mb-4"),
                    html.P("Watch your bacterial colony evolve in real-time",
                           className="text-center text-muted mb-4")
                ])
            ]),
            
            # Control Panel
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("Simulation Controls"),
                        dbc.CardBody([
                            dbc.ButtonGroup([
                                dbc.Button("▶ Play", 
                                         id="play-button", 
                                         color="success",
                                         className="me-2"),
                                dbc.Button("⏸ Pause",
                                         id="pause-button",
                                         color="warning",
                                         className="me-2"),
                                dbc.Button("↺ Reset", 
                                         id="reset-button", 
                                         color="danger")
                            ], className="mb-3"),
                            html.Label("Simulation Speed"),
                            dcc.Slider(
                                id='speed-slider',
                                min=0.1,
                                max=2.0,
                                step=0.1,
                                value=1.0,
                                marks={i/10: str(i/10) for i in range(1, 21)},
                                className="mb-3"
                            ),
                            html.Div(id='current-step-display',
                                    className="text-muted")
                        ])
                    ], className="mb-4")
                ])
            ]),
            
            # Main content
            dbc.Row([
                # Left column: Simulation view and stats
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("Simulation View"),
                        dbc.CardBody([
                            dcc.Graph(id='simulation-graph')
                        ])
                    ], className="mb-4"),
                    
                    dbc.Card([
                        dbc.CardHeader("Current Statistics"),
                        dbc.CardBody([
                            html.Div(id='statistics-display')
                        ])
                    ])
                ], width=6),
                
                # Right column: Time series
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("Time Series Data"),
                        dbc.CardBody([
                            dcc.Graph(id='time-series-graph')
                        ])
                    ])
                ], width=6)
            ]),
            
            # Hidden components
            dcc.Store(id='simulation-state', data=self.initial_state),
            dcc.Interval(
                id='interval-component',
                interval=200,  # Default 200ms interval
                n_intervals=0,
                disabled=True
            )
        ], fluid=True, className="p-4")

    def setup_callbacks(self):
        """Setup the dashboard callbacks"""
        @self.app.callback(
            Output('interval-component', 'disabled'),
            [Input('play-button', 'n_clicks'),
             Input('pause-button', 'n_clicks'),
             Input('reset-button', 'n_clicks')],
            [State('interval-component', 'disabled')]
        )
        def control_simulation(play, pause, reset, current_disabled):
            ctx = dash.callback_context
            if not ctx.triggered:
                return True
            
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            
            if button_id == 'play-button':
                return False
            elif button_id == 'pause-button':
                return True
            elif button_id == 'reset-button':
                return True
            return current_disabled

        @self.app.callback(
            [Output('simulation-graph', 'figure'),
             Output('time-series-graph', 'figure'),
             Output('statistics-display', 'children'),
             Output('current-step-display', 'children'),
             Output('simulation-state', 'data')],
            [Input('interval-component', 'n_intervals'),
             Input('reset-button', 'n_clicks')],
            [State('simulation-state', 'data')]
        )
        def update_graphs(n_intervals, reset_clicks, current_state):
            ctx = dash.callback_context
            triggered_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None
            
            if triggered_id == 'reset-button' and reset_clicks is not None and reset_clicks > self.reset_clicks:
                self.current_step = 0
                self.reset_histories()
                self.reset_clicks = reset_clicks
            elif triggered_id == 'interval-component' and n_intervals is not None:
                self.current_step = min(self.current_step + 1, self.max_steps - 1)
            
            # Get current state
            ecoli_list, glucose_list, metabolite_list = self.states[self.current_step]
            
            # Update histories
            self.time_history.append(self.current_step)
            self.metabolites_history.append(len(metabolite_list))
            self.population_history.append(len(ecoli_list))
            self.atp_history.append(np.mean([e.ATP for e in ecoli_list]) if ecoli_list else 0)
            self.excreted_metabolites_history.append(sum(e.produced_metabolite for e in ecoli_list))
            self.toxic_concentration_history.append(self.calculate_average_toxicity(ecoli_list, metabolite_list))
            
            # Create views
            sim_fig = self.create_simulation_view(ecoli_list, glucose_list, metabolite_list)
            ts_fig = self.create_time_series_view()
            stats = self.create_statistics_display(ecoli_list, glucose_list, metabolite_list)
            step_display = f"Current Step: {self.current_step + 1}/{self.max_steps}"
            
            # Update state
            current_state = {
                'current_step': self.current_step,
                'histories': {
                    'time': list(self.time_history),
                    'metabolites': list(self.metabolites_history),
                    'population': list(self.population_history),
                    'atp': list(self.atp_history),
                    'excreted_metabolites': list(self.excreted_metabolites_history),
                    'toxic_concentration': list(self.toxic_concentration_history)
                }
            }
            
            return sim_fig, ts_fig, stats, step_display, current_state
        
        @self.app.callback(
            Output('interval-component', 'interval'),
            [Input('speed-slider', 'value')]
        )
        def update_interval(value):
            return int(200 / value)  # Base interval is 200ms for smoother animation

    def create_simulation_view(self, ecoli_list, glucose_list, metabolite_list):
        """Create the main simulation view with all particles"""
        fig = go.Figure()
        
        # Add drain visualization
        for start_x, end_x in self.drain_system.drain_positions:
            # Add semi-transparent rectangle for drain
            fig.add_shape(
                type="rect",
                x0=start_x * self.space_size * 1e6,
                x1=end_x * self.space_size * 1e6,
                y0=0,
                y1=self.drain_system.drain_height * self.space_size * 1e6,
                fillcolor="rgba(128, 128, 128, 0.3)",
                line=dict(color="rgba(128, 128, 128, 0.5)"),
                layer="below"
            )
        
        # Add metabolites with concentration-based visualization
        if metabolite_list:
            fig.add_trace(go.Scatter(
                x=[m.x * 1e6 for m in metabolite_list],
                y=[m.y * 1e6 for m in metabolite_list],
                mode='markers',
                name='Metabolites',
                marker=dict(
                    size=4,
                    color=[m.concentration for m in metabolite_list],
                    colorscale='RdPu',  # Red-Purple scale for toxicity
                    opacity=0.4,
                    colorbar=dict(
                        title='Metabolite<br>Concentration',
                        thickness=20,  # Increased thickness
                        len=0.6,  # Increased length
                        titlefont=dict(size=14),  # Larger title font
                        tickfont=dict(size=12)  # Larger tick font
                    )
                )
            ))
        
        # Add glucose particles
        fig.add_trace(go.Scatter(
            x=[g.x * 1e6 for g in glucose_list],
            y=[g.y * 1e6 for g in glucose_list],
            mode='markers',
            name='Glucose',
            marker=dict(
                size=5,
                color='#17a2b8',
                opacity=0.6
            )
        ))
        
        # Add bacteria with ATP-based coloring
        fig.add_trace(go.Scatter(
            x=[e.x * 1e6 for e in ecoli_list],
            y=[e.y * 1e6 for e in ecoli_list],
            mode='markers',
            name='E. coli',
            marker=dict(
                size=8,
                color=[min(e.ATP/200, 1) for e in ecoli_list],
                colorscale='Viridis',
                opacity=0.8,
                colorbar=dict(
                    title='Bacterial<br>ATP',
                    thickness=20,  # Increased thickness
                    len=0.6,  # Increased length
                    x=1.15,  # Moved slightly right
                    titlefont=dict(size=14),  # Larger title font
                    tickfont=dict(size=12)  # Larger tick font
                )
            )
        ))
        
        # Update layout with larger legend
        fig.update_layout(
            template="plotly_dark",
            margin=dict(l=40, r=140, t=40, b=40),  # Increased right margin for colorbar
            xaxis=dict(
                title='X Position (µm)',
                showgrid=True,
                gridcolor='rgba(255, 255, 255, 0.1)',
                titlefont=dict(size=14),  # Larger axis title
                tickfont=dict(size=12)  # Larger tick labels
            ),
            yaxis=dict(
                title='Y Position (µm)',
                showgrid=True,
                gridcolor='rgba(255, 255, 255, 0.1)',
                titlefont=dict(size=14),  # Larger axis title
                tickfont=dict(size=12)  # Larger tick labels
            ),
            showlegend=True,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.85,  # Adjusted for colorbar
                font=dict(size=14),  # Larger legend text
                bgcolor="rgba(0,0,0,0.3)",  # Semi-transparent background
                bordercolor="rgba(255,255,255,0.2)",
                borderwidth=1
            )
        )
        
        return fig

    def create_time_series_view(self):
        """Create time series plots for various metrics"""
        fig = make_subplots(
            rows=5, cols=1,
            subplot_titles=(
                'Environmental Metabolites Over Time',
                'Bacterial Population Over Time',
                'Average ATP Levels Over Time',
                'Total Metabolites Produced Over Time',
                'Average Toxic Concentration Over Time'
            ),
            vertical_spacing=0.12  # Increased spacing between subplots
        )
        
        # Add time series traces with custom styling and clear names
        fig.add_trace(
            go.Scatter(
                x=list(self.time_history), 
                y=list(self.metabolites_history),
                name='Active Metabolites',
                line=dict(color='#28a745', width=2),
                showlegend=True
            ),
            row=1, col=1
        )
        
        fig.add_trace(
            go.Scatter(
                x=list(self.time_history), 
                y=list(self.population_history),
                name='Bacterial Count',
                line=dict(color='#dc3545', width=2),
                showlegend=True
            ),
            row=2, col=1
        )
        
        fig.add_trace(
            go.Scatter(
                x=list(self.time_history), 
                y=list(self.atp_history),
                name='Average ATP',
                line=dict(color='#ffc107', width=2),
                showlegend=True
            ),
            row=3, col=1
        )
        
        fig.add_trace(
            go.Scatter(
                x=list(self.time_history), 
                y=list(self.excreted_metabolites_history),
                name='Total Metabolites',
                line=dict(color='#17a2b8', width=2),
                showlegend=True
            ),
            row=4, col=1
        )
        
        fig.add_trace(
            go.Scatter(
                x=list(self.time_history), 
                y=list(self.toxic_concentration_history),
                name='Toxicity Level',
                line=dict(color='#e83e8c', width=2),
                showlegend=True
            ),
            row=5, col=1
        )
        
        # Update layout with improved styling
        fig.update_layout(
            height=1200,  # Set height here instead
            showlegend=True,
            template="plotly_dark",
            margin=dict(l=50, r=50, t=60, b=50),
            legend=dict(
                orientation="h",  # Horizontal legend
                yanchor="bottom",
                y=-0.2,  # Position below the plots
                xanchor="center",
                x=0.5,
                bgcolor="rgba(0,0,0,0.3)",  # Semi-transparent background
                bordercolor="rgba(255,255,255,0.2)",
                borderwidth=1
            )
        )
        
        # Update axes with improved labels
        for i in range(1, 6):
            fig.update_xaxes(
                title_text="Time Steps",
                showgrid=True,
                gridwidth=1,
                gridcolor='rgba(255, 255, 255, 0.1)',
                row=i, col=1
            )
        
        # Add specific y-axis labels
        y_labels = [
            "Count",
            "Count",
            "ATP Units",
            "Count",
            "Concentration"
        ]
        
        for i, label in enumerate(y_labels, 1):
            fig.update_yaxes(
                title_text=label,
                showgrid=True,
                gridwidth=1,
                gridcolor='rgba(255, 255, 255, 0.1)',
                row=i, col=1
            )
        
        return fig

    def create_statistics_display(self, ecoli_list, glucose_list, metabolite_list):
        """Create the statistics display panel"""
        # Calculate average toxic concentration
        avg_toxic_conc = self.calculate_average_toxicity(ecoli_list, metabolite_list)
        
        stats = [
            dbc.Row([
                dbc.Col([
                    html.H3(str(len(ecoli_list)), className="text-primary"),
                    html.P("Bacteria", className="text-muted")
                ], className="text-center mb-3"),
                
                dbc.Col([
                    html.H3(str(len(glucose_list)), className="text-info"),
                    html.P("Glucose Particles", className="text-muted")
                ], className="text-center mb-3")
            ]),
            dbc.Row([
                dbc.Col([
                    html.H3(str(len(metabolite_list)), className="text-success"),
                    html.P("Active Metabolites", className="text-muted")
                ], className="text-center mb-3"),
                
                dbc.Col([
                    html.H3(f"{np.mean([e.ATP for e in ecoli_list]):.2f}" if ecoli_list else "0", 
                           className="text-warning"),
                    html.P("Average ATP", className="text-muted")
                ], className="text-center")
            ]),
            dbc.Row([
                dbc.Col([
                    html.H3(f"{sum(e.produced_metabolite for e in ecoli_list)}", 
                           className="text-info"),
                    html.P("Total Metabolites Produced", className="text-muted")
                ], className="text-center mb-3"),
                
                dbc.Col([
                    html.H3(f"{avg_toxic_conc:.2f}",
                           className="text-danger"),
                    html.P("Avg. Toxic Concentration", className="text-muted")
                ], className="text-center")
            ]),
            dbc.Row([
                dbc.Col([
                    html.H3(f"{self.drain_system.collected_metabolites}",
                           className="text-secondary"),
                    html.P("Collected Metabolites", className="text-muted")
                ], className="text-center")
            ])
        ]
        return stats

# ====================
# --- Main ---
# ====================

def run_multiple_simulations(scenarios):
    """
    Run multiple simulations with different parameters and compare results
    
    Args:
        scenarios (list): List of dictionaries containing simulation parameters
        
    Returns:
        list: Results from each simulation
    """
    results = []
    
    for scenario in scenarios:
        print(f"\nRunning scenario: {scenario['name']}")
        print("Parameters:", {k: v for k, v in scenario.items() if k != 'name'})
        
        # Run simulation with scenario parameters
        states, space_size, drain_system = simulate(
            n_ecoli=scenario.get('n_ecoli', 50),
            n_glucose=scenario.get('n_glucose', 1000),
            space_size=scenario.get('space_size', 20e-6),
            max_steps=scenario.get('max_steps', 250),
            glucose_feed_interval=scenario.get('glucose_feed_interval', 20),
            glucose_feed_amount=scenario.get('glucose_feed_amount', 50)
        )
        
        # Calculate final statistics
        final_ecoli = states[-1][0]
        final_glucose = states[-1][1]
        final_metabolites = states[-1][2]
        
        # Collect results
        result = {
            'scenario': scenario['name'],
            'final_bacteria': len(final_ecoli),
            'final_glucose': len(final_glucose),
            'final_metabolites': len(final_metabolites),
            'collected_metabolites': drain_system.collected_metabolites,
            'avg_final_atp': np.mean([e.ATP for e in final_ecoli]) if final_ecoli else 0,
            'total_steps': len(states)
        }
        results.append(result)
        
        print(f"\nResults for {scenario['name']}:")
        for key, value in result.items():
            if key != 'scenario':
                print(f"{key}: {value}")
    
    return results

if __name__ == '__main__':
    # Define available scenarios
    scenarios = {
        '1': {
            'name': 'baseline',
            'description': 'Standard simulation parameters',
            'n_ecoli': 50,
            'n_glucose': 1000,
            'max_steps': 500
        },
        '2': {
            'name': 'high_density',
            'description': 'Higher bacterial and glucose density in a smaller space',
            'n_ecoli': 100,
            'n_glucose': 2000,
            'space_size': 40e-6,
            'max_steps': 500
        },
        '3': {
            'name': 'low_nutrient',
            'description': 'Limited nutrients with slower feeding rate',
            'n_ecoli': 50,
            'n_glucose': 500,
            'glucose_feed_interval': 40,
            'glucose_feed_amount': 25,
            'max_steps': 500
        },
        '4': {
            'name': 'efficient_drainage',
            'description': 'Optimized for metabolite collection',
            'n_ecoli': 75,
            'n_glucose': 1500,
            'max_steps': 500
        }
    }
    
    # Print available scenarios
    print("\nAvailable Scenarios:")
    print("-" * 50)
    for key, scenario in scenarios.items():
        print(f"{key}. {scenario['name']}")
        print(f"   {scenario['description']}")
        print(f"   Parameters: {', '.join(f'{k}: {v}' for k, v in scenario.items() if k not in ['name', 'description'])}")
        print()
    
    # Get user choice
    while True:
        choice = input("\nSelect a scenario (1-4): ").strip()
        if choice in scenarios:
            break
        print("Invalid choice. Please select a number between 1 and 4.")
    
    # Get selected scenario
    selected_scenario = scenarios[choice]
    print(f"\nRunning scenario: {selected_scenario['name']}")
    print("Parameters:", {k: v for k, v in selected_scenario.items() if k not in ['name', 'description']})
    
    # Get space size from scenario or use default
    space_size = selected_scenario.get('space_size', 20e-6)
    
    # Initialize stream field for fluid flow
    fluid = StreamField(width=space_size, height=space_size)
    
    # Run simulation
    states, space_size, drain_system = simulate(
        n_ecoli=selected_scenario['n_ecoli'],
        n_glucose=selected_scenario['n_glucose'],
        space_size=space_size,
        max_steps=selected_scenario.get('max_steps', 250),
        glucose_feed_interval=selected_scenario.get('glucose_feed_interval', 20),
        glucose_feed_amount=selected_scenario.get('glucose_feed_amount', 50)
    )
    
    # Export simulation data to CSV
    csv_filename = export_detailed_data(states)
    print(f"\nSimulation data has been exported to: {csv_filename}")

    #export simulation data to PNG
    export_summary_image(states)
    
    # Initialize and run dashboard
    dashboard = SimulationDashboard(states, space_size, fluid_field=fluid, drain_system=drain_system)
    dashboard.run(debug=False, port=8050)
