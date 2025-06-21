import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from collections import deque
import copy
import warnings
import os
from datetime import datetime
import csv

# Suppress deprecation warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ==========================
# --- Configuration ---
# ==========================

@dataclass
class SimulationConfig:
    """A central class to hold all simulation parameters."""
    # General Simulation
    SPACE_SIZE: float = 20e-6

    # E. coli Parameters
    ECOLI_INITIAL_ATP: float = 200.0
    ECOLI_WIDTH: float = 2e-6
    ECOLI_HEIGHT: float = 1e-6
    ECOLI_GRAVITY: float = 4.5e-6
    ECOLI_GLUCOSE_SENSITIVITY: float = 0.5
    ECOLI_BASE_GROWTH_RATE: float = 0.12
    ECOLI_ACETATE_CONSUMPTION_RADIUS: float = 2e-6
    ECOLI_SENSING_RADIUS: float = 5e-6
    ECOLI_CHEMOTAXIS_STEP_SIZE: float = 1e-6
    ECOLI_CHEMOTAXIS_NOISE: float = 0.5
    ECOLI_CHEMOTAXIS_DETECTION_RADIUS: float = 1e-5
    ECOLI_DIVISION_THRESHOLD_ATP: float = 200.0
    ECOLI_SURVIVAL_THRESHOLD_ATP: float = 10.0
    ECOLI_METABOLISM_ATP_COST: float = 0.5
    ECOLI_DEATH_PROB_BASE: float = 0.005
    ECOLI_DEATH_PROB_LOW_ATP: float = 0.1
    ECOLI_DEATH_PROB_VERY_LOW_ATP: float = 0.2
    ECOLI_DEATH_PROB_CRITICAL_ATP: float = 0.5
    ECOLI_LOW_ATP_THRESHOLD: float = 50.0
    ECOLI_VERY_LOW_ATP_THRESHOLD: float = 25.0
    ECOLI_CRITICAL_ATP_THRESHOLD: float = 10.0
    ECOLI_FLUID_H_FACTOR: float = 15.0
    ECOLI_FLUID_V_FACTOR: float = 10.0
    ECOLI_ATP_GAIN_PER_GLUCOSE: float = 10.0

    # Glucose Parameters
    GLUCOSE_DIAMETER: float = 0.0009e-6
    GLUCOSE_GRAVITY: float = 0.0025e-6
    GLUCOSE_DIFFUSION_BASE: float = 2.5e-6
    GLUCOSE_MOMENTUM_DECAY: float = 0.95
    GLUCOSE_BOUNCE_FACTOR: float = 0.7
    GLUCOSE_FLUID_INTERACTION_FACTOR: float = 2.0

    # Acetate (Metabolite) Parameters
    ACETATE_DIAMETER: float = 0.0005e-6
    ACETATE_DIFFUSION_RATE: float = 2.5e-7
    ACETATE_DECAY_RATE: float = 0.0
    ACETATE_INITIAL_CONCENTRATION: float = 0.3
    ACETATE_TOXICITY_THRESHOLD: float = 8.0
    ACETATE_LETHAL_CONCENTRATION: float = 60.0
    ACETATE_MAX_TOXICITY_EFFECT: float = 1.0
    ACETATE_GRAVITY: float = 0.01e-6
    ACETATE_DEATH_PROB_FACTOR: float = 0.8
    ACETATE_ATP_DRAIN_FACTOR: float = 0.4

    # Fluid Dynamics Parameters
    FLUID_BASE_DOWNWARD_FLOW: float = -0.1e-6
    FLUID_DRAIN_STRENGTH: float = -0.2e-6
    FLUID_UPWARD_FLOW_STRENGTH: float = 0.2e-6
    FLUID_DRAIN_POSITIONS: list = field(default_factory=lambda: [(0.3, 0.4), (0.6, 0.7)])
    
    # Drain System Parameters
    DRAIN_HEIGHT_PERCENTAGE: float = 0.08
    DRAIN_COLLECTION_RADIUS_PERCENTAGE: float = 0.02
    DRAIN_RECYCLING_DELAY_STEPS: int = 25

    # Spatial Grid Parameters
    BACTERIA_GRID_CELL_SIZE: float = 5e-6
    GLUCOSE_GRID_CELL_SIZE: float = 2e-6
    METABOLITE_GRID_CELL_SIZE: float = 3e-6


# Dash imports with error handling
try:
    import dash
    from dash import dcc, html
    from dash.dependencies import Input, Output, State
    import plotly.graph_objs as go
    from plotly.subplots import make_subplots
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
    def __init__(self, width, height, config):
        self.width = width
        self.height = height
        self.config = config
        
        # Adjusted flow parameters for better suspension
        self.base_flow = config.FLUID_BASE_DOWNWARD_FLOW
        self.drain_strength = config.FLUID_DRAIN_STRENGTH
        self.upward_flow_strength = config.FLUID_UPWARD_FLOW_STRENGTH
        self.drain_positions = config.FLUID_DRAIN_POSITIONS
        
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

def create_spatial_grids(config):
    """Create spatial grids for different particle types"""
    space_size = config.SPACE_SIZE
    # Use different cell sizes based on typical search radii
    bacteria_grid = SpatialGrid(space_size, space_size, cell_size=config.BACTERIA_GRID_CELL_SIZE)
    glucose_grid = SpatialGrid(space_size, space_size, cell_size=config.GLUCOSE_GRID_CELL_SIZE)
    metabolite_grid = SpatialGrid(space_size, space_size, cell_size=config.METABOLITE_GRID_CELL_SIZE)
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
    config: SimulationConfig
    collected_metabolites: int = 0
    drain_positions: list = field(default_factory=lambda: [(0.3, 0.4), (0.6, 0.7)])
    drain_height: float = field(default=0.08)
    collection_radius: float = field(default=0.02)
    recycling_queue: list = field(default_factory=list)
    recycling_delay: int = field(default=25)
    
    def __post_init__(self):
        self.drain_positions = self.config.FLUID_DRAIN_POSITIONS
        self.drain_height = self.config.DRAIN_HEIGHT_PERCENTAGE
        self.collection_radius = self.config.DRAIN_COLLECTION_RADIUS_PERCENTAGE
        self.recycling_delay = self.config.DRAIN_RECYCLING_DELAY_STEPS
    
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
    config: SimulationConfig
    type: str = "generic"
    vertical_momentum: float = 0.0
    horizontal_momentum: float = 0.0
    diameter: float = field(init=False)
    diffusion_rate: float = field(init=False)
    decay_rate: float = field(init=False)
    concentration: float = field(init=False)
    toxicity_threshold: float = field(init=False)
    lethal_concentration: float = field(init=False)
    max_toxicity: float = field(init=False)
    gravity: float = field(init=False)

    def __post_init__(self):
        self.diameter = self.config.ACETATE_DIAMETER
        self.diffusion_rate = self.config.ACETATE_DIFFUSION_RATE
        self.decay_rate = self.config.ACETATE_DECAY_RATE
        self.concentration = self.config.ACETATE_INITIAL_CONCENTRATION
        self.toxicity_threshold = self.config.ACETATE_TOXICITY_THRESHOLD
        self.lethal_concentration = self.config.ACETATE_LETHAL_CONCENTRATION
        self.max_toxicity = self.config.ACETATE_MAX_TOXICITY_EFFECT
        self.gravity = self.config.ACETATE_GRAVITY
    
    def calculate_toxicity(self, local_concentration):
        """
        Calculate toxic effect based on local concentration using a linear model.
        Toxicity starts at `toxicity_threshold` and reaches `max_toxicity` at `lethal_concentration`.
        """
        if local_concentration <= self.toxicity_threshold:
            return 0.0
        
        if local_concentration >= self.lethal_concentration:
            return self.max_toxicity

        # Linear interpolation
        toxic_range = self.lethal_concentration - self.toxicity_threshold
        if toxic_range <= 0:
            return self.max_toxicity  # Avoid division by zero
            
        toxicity_level = (local_concentration - self.toxicity_threshold) / toxic_range
        return self.max_toxicity * toxicity_level
    
    def diffuse(self, diffusion_strength, space_size, fluid_field=None):
        """Diffuse metabolite using a momentum-based system similar to glucose."""
        momentum_decay = 0.98  # Slightly less decay than glucose
        bounce_factor = 0.5    # Less bouncy than glucose
        
        # Update momentum with random forces
        self.horizontal_momentum = (self.horizontal_momentum * momentum_decay + 
                                  np.random.normal(0, self.diffusion_rate))
        self.vertical_momentum = (self.vertical_momentum * momentum_decay + 
                                np.random.normal(0, self.diffusion_rate))
        
        # Apply momentum to position
        self.x += self.horizontal_momentum * diffusion_strength
        self.y += self.vertical_momentum * diffusion_strength
        
        # Apply gravity (stronger effect)
        self.y -= self.gravity * diffusion_strength
        
        # Apply fluid effects to momentum
        if fluid_field is not None:
            u, v = fluid_field.get_velocity_at_point(self.x, self.y)
            fluid_factor = self.config.GLUCOSE_FLUID_INTERACTION_FACTOR * 0.5 # Less sensitive than glucose
            
            self.horizontal_momentum += u * diffusion_strength * fluid_factor
            self.vertical_momentum += v * diffusion_strength * fluid_factor
        
        # Boundary interactions
        if self.x <= 0:
            self.x = 0
            self.horizontal_momentum = abs(self.horizontal_momentum) * bounce_factor
        elif self.x >= space_size:
            self.x = space_size
            self.horizontal_momentum = -abs(self.horizontal_momentum) * bounce_factor
            
        if self.y <= 0:
            self.y = 0
            self.vertical_momentum = abs(self.vertical_momentum) * bounce_factor * 0.1 # Stick to bottom
        elif self.y >= space_size:
            self.y = space_size
            self.vertical_momentum = -abs(self.vertical_momentum) * bounce_factor
        
        # Increase concentration decay over time
        self.concentration *= (1 - self.decay_rate * 1.2)

    def __deepcopy__(self, memo):
        new_metabolite = Metabolite(
            x=self.x,
            y=self.y,
            config=self.config,
            type=self.type,
        )
        new_metabolite.concentration = self.concentration
        new_metabolite.vertical_momentum = self.vertical_momentum
        new_metabolite.horizontal_momentum = self.horizontal_momentum
        return new_metabolite

@dataclass
class Glucose:
    x: float
    y: float
    config: SimulationConfig
    diameter: float = field(init=False)
    gravity: float = field(init=False)
    vertical_momentum: float = 0.0
    horizontal_momentum: float = 0.0

    def __post_init__(self):
        self.diameter = self.config.GLUCOSE_DIAMETER
        self.gravity = self.config.GLUCOSE_GRAVITY

    def diffuse(self, diffusion_strength, space_size, fluid_field=None):
        # Keep existing diffusion parameters but increase gravity effect
        base_diffusion = self.config.GLUCOSE_DIFFUSION_BASE
        momentum_decay = self.config.GLUCOSE_MOMENTUM_DECAY
        
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
            fluid_factor = self.config.GLUCOSE_FLUID_INTERACTION_FACTOR
            
            self.horizontal_momentum += u * diffusion_strength * fluid_factor
            self.vertical_momentum += v * diffusion_strength * fluid_factor * 0.3
        
        # Boundary interactions
        bounce_factor = self.config.GLUCOSE_BOUNCE_FACTOR
        
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
        new_glucose = Glucose(x=self.x, y=self.y, config=self.config)
        new_glucose.vertical_momentum = self.vertical_momentum
        new_glucose.horizontal_momentum = self.horizontal_momentum
        return new_glucose

@dataclass
class Ecoli:
    x: float
    y: float
    config: SimulationConfig
    ATP: float = field(init=False)
    width: float = field(init=False)
    height: float = field(init=False)
    produced_metabolite: int = 0
    gravity: float = field(init=False)
    glucose_sensitivity: float = field(init=False)
    base_growth_rate: float = field(init=False)
    metabolite_sensing_radius: float = field(init=False)
    
    def __post_init__(self):
        self.ATP = self.config.ECOLI_INITIAL_ATP
        self.width = self.config.ECOLI_WIDTH
        self.height = self.config.ECOLI_HEIGHT
        self.gravity = self.config.ECOLI_GRAVITY
        self.glucose_sensitivity = self.config.ECOLI_GLUCOSE_SENSITIVITY
        self.base_growth_rate = self.config.ECOLI_BASE_GROWTH_RATE
        self.metabolite_sensing_radius = self.config.ECOLI_SENSING_RADIUS
    
    def calculate_growth_probability(self, glucose_list, glucose_grid):
        """Calculate growth probability using spatial grid for efficiency"""
        nearby_indices = glucose_grid.get_nearby_indices(self.x, self.y, self.metabolite_sensing_radius)
        # Filter out invalid indices
        valid_indices = [i for i in nearby_indices if i < len(glucose_list)]
        local_glucose = [glucose_list[i] for i in valid_indices if glucose_list[i] is not None]
        
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
                local_concentration += metabolite.concentration * weight
        
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
            atp_reduction = self.ATP * max_toxic_effect * self.config.ACETATE_ATP_DRAIN_FACTOR
            self.ATP -= atp_reduction
            return atp_reduction
        return 0

    def excrete_metabolite(self, space_size, glucose_list, glucose_grid):
        """Excrete acetate if enough glucose is available to be consumed."""
        consumption_radius = self.config.ECOLI_ACETATE_CONSUMPTION_RADIUS
        # Find nearby glucose that has not been consumed yet
        nearby_indices = glucose_grid.get_nearby_indices(self.x, self.y, consumption_radius)
        
        available_glucose = []
        for i in nearby_indices:
            if i < len(glucose_list) and glucose_list[i] is not None:
                g = glucose_list[i]
                distance = np.hypot(self.x - g.x, self.y - g.y)
                if distance < consumption_radius:
                    available_glucose.append((i, g))

        # If there are at least 2 glucose molecules nearby, consume them to produce acetate
        if len(available_glucose) >= 2:
            # Sort by distance to consume the closest ones
            available_glucose.sort(key=lambda item: np.hypot(self.x - item[1].x, self.y - item[1].y))
            
            indices_to_consume = [item[0] for item in available_glucose[:2]]
            
            # Create new metabolite (acetate) near the bacteria
            offset_x = np.random.uniform(-self.width/2, self.width/2)
            offset_y = np.random.uniform(-self.height/2, self.height/2)
            
            metabolite = Metabolite(
                x=np.clip(self.x + offset_x, 0, space_size),
                y=np.clip(self.y + offset_y, 0, space_size),
                config=self.config,
                type="acetate"
            )
            
            self.produced_metabolite += 1
            
            # Increase ATP when consuming glucose
            self.ATP += self.config.ECOLI_ATP_GAIN_PER_GLUCOSE * 2
            
            return metabolite, indices_to_consume
            
        return None, []

    def move(self, glucose_list, glucose_grid, space_size, fluid_field=None):
        """Stochastic movement with chemotaxis"""
        # Get nearby glucose for chemotaxis
        _, nearby_glucose = self.calculate_growth_probability(glucose_list, glucose_grid)
        
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
        
        step_size = self.config.ECOLI_CHEMOTAXIS_STEP_SIZE
        noise = self.config.ECOLI_CHEMOTAXIS_NOISE

        # Combine directed movement, noise, and very gentle gravity
        final_dx = step_size * ((1 - noise) * dx + noise * ndx)
        final_dy = step_size * ((1 - noise) * dy + noise * ndy) - self.gravity * step_size * 2  # Reduced gravity effect
        
        # Apply fluid velocity if fluid field exists
        if fluid_field is not None:
            u, v = fluid_field.get_velocity_at_point(self.x, self.y)
            # Use configurable factors for fluid interaction
            final_dx += u * step_size * self.config.ECOLI_FLUID_H_FACTOR
            final_dy += v * step_size * self.config.ECOLI_FLUID_V_FACTOR
        
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

    def can_divide(self):
        """Check if cell has enough ATP to divide"""
        return self.ATP >= self.config.ECOLI_DIVISION_THRESHOLD_ATP

    def divide(self, space_size):
        """Discrete cell division event"""
        self.ATP /= 2
        # New bacteria appears slightly to the side and above/below
        offset_x = np.random.uniform(-1e-6, 1e-6)
        offset_y = np.random.uniform(-1e-6, 1e-6)
        return Ecoli(
            x=np.clip(self.x + offset_x, 0, space_size),
            y=np.clip(self.y + offset_y, 0, space_size),
            config=self.config
        )

    def is_alive(self):
        """Check if cell has enough ATP to stay alive"""
        return self.ATP >= self.config.ECOLI_SURVIVAL_THRESHOLD_ATP

    def __deepcopy__(self, memo):
        new_ecoli = Ecoli(
            x=self.x,
            y=self.y,
            config=self.config
        )
        new_ecoli.ATP = self.ATP
        new_ecoli.produced_metabolite = self.produced_metabolite
        return new_ecoli

# ====================
# --- Gillespie algorithm ---
# ====================

def gillespie_step(ecoli, glucose_list, glucose_grid, metabolite_list, metabolite_grid, space_size):
    """
    Gillespie algorithm step for discrete event simulation.
    Events: division, death, death by toxicity, or idle.
    """
    config = ecoli.config
    # Calculate growth probability based on local glucose
    growth_probability, _ = ecoli.calculate_growth_probability(glucose_list, glucose_grid)
    
    # ATP decay over time (bacteria consume energy to stay alive)
    ecoli.ATP -= config.ECOLI_METABOLISM_ATP_COST
    
    # Calculate death probability based on ATP level
    death_probability = config.ECOLI_DEATH_PROB_BASE
    if ecoli.ATP < config.ECOLI_LOW_ATP_THRESHOLD:
        death_probability = config.ECOLI_DEATH_PROB_LOW_ATP
        if ecoli.ATP < config.ECOLI_VERY_LOW_ATP_THRESHOLD:
            death_probability = config.ECOLI_DEATH_PROB_VERY_LOW_ATP
        if ecoli.ATP < config.ECOLI_CRITICAL_ATP_THRESHOLD:
            death_probability = config.ECOLI_DEATH_PROB_CRITICAL_ATP
            
    # Calculate toxicity death probability
    local_concentration = ecoli.calculate_local_metabolite_concentration(metabolite_list, metabolite_grid)
    dummy_metabolite = Metabolite(x=0, y=0, config=ecoli.config)  # To access calculate_toxicity
    toxic_effect = dummy_metabolite.calculate_toxicity(local_concentration)
    toxicity_death_prob = toxic_effect * config.ACETATE_DEATH_PROB_FACTOR
    
    # Define possible events and their rates
    rates = {
        'divide': growth_probability if ecoli.can_divide() else 0,
        'die': death_probability,
        'die_by_toxicity': toxicity_death_prob,
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

def add_glucose_from_top(glucose_list, n_glucose, config):
    """
    Add new glucose molecules with improved floating behavior
    """
    space_size = config.SPACE_SIZE
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
            glucose = Glucose(x_position, y_position, config=config)
            glucose.vertical_momentum = np.random.uniform(0.5e-6, 1.5e-6)  # Initial upward momentum
            glucose.horizontal_momentum = np.random.normal(0, 1e-6)  # Random horizontal momentum
            
            new_glucose.append(glucose)
    
    glucose_list.extend(new_glucose)
    return glucose_list

def create_ecoli_and_glucose(n_ecoli, n_glucose, config):
    space_size = config.SPACE_SIZE
    ecoli_list = [Ecoli(np.random.uniform(0, space_size), np.random.uniform(0, space_size), config)
                  for _ in range(n_ecoli)]
    glucose_list = [Glucose(np.random.uniform(0, space_size), np.random.uniform(0, space_size), config)
                    for _ in range(n_glucose)]
    return ecoli_list, glucose_list

def simulate(n_ecoli=10, n_glucose=300, space_size=20e-6, max_steps=250, 
            glucose_feed_interval=20, glucose_feed_amount=50, use_drain=True, config=None):
    """
    Simulate bacterial growth with simplified fluid dynamics and discrete events
    """
    # Initialize config if not provided
    if config is None:
        config = SimulationConfig(SPACE_SIZE=space_size)

    # Initialize spatial grids
    bacteria_grid, glucose_grid, metabolite_grid = create_spatial_grids(config)
    
    # Initialize stream field for fluid flow only if the drain is active
    fluid = StreamField(width=space_size, height=space_size, config=config) if use_drain else None
    
    # Initialize drain system
    drain_system = DrainSystem(space_size=space_size, config=config) if use_drain else None
    
    ecoli_list, glucose_list = create_ecoli_and_glucose(n_ecoli, n_glucose, config)
    metabolite_list = []
    states = []
    total_metabolites_produced = 0
    toxicity_deaths = 0

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
        collected_count = drain_system.collected_metabolites if drain_system else 0
        print(f"\rProgress: {progress:3.1f}% | Step: {step}/{max_steps} | "
              f"Bacteria: {len(ecoli_list)} | Glucose: {len(glucose_list)} | "
              f"Metabolites: {len(metabolite_list)} | Collected: {collected_count}", 
              end="", flush=True)
        
        # Process recycling
        if drain_system:
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
                if not (drain_system and drain_system.check_collection(metabolite, "metabolite")):
                    # Apply diffusion before adding to new list
                    metabolite.diffuse(diffusion_strength=1.0, space_size=space_size, fluid_field=fluid)
                    new_metabolites.append(metabolite)
            metabolite_list = new_metabolites
        
        # Process glucose with drain collection
        if glucose_list:
            new_glucose_list = []
            for glucose in glucose_list:
                if not (drain_system and drain_system.check_collection(glucose, "glucose")):
                    # Apply diffusion before adding to new list
                    glucose.diffuse(diffusion_strength=1.0, space_size=space_size, fluid_field=fluid)
                    new_glucose_list.append(glucose)
            glucose_list = new_glucose_list
        
        # Process bacteria with drain collection
        new_bact = []
        survivors = []
        for bacteria in ecoli_list:
            if not (drain_system and drain_system.check_collection(bacteria, "bacteria")):
                # Move bacteria
                bacteria.move(glucose_list, glucose_grid, space_size, fluid_field=fluid)

                # Apply non-lethal toxic effects (ATP drain)
                bacteria.apply_toxic_effects(metabolite_list, metabolite_grid)
                
                # Handle metabolite excretion
                new_metabolite, consumed_indices = bacteria.excrete_metabolite(space_size, glucose_list, glucose_grid)
                if new_metabolite:
                    metabolite_list.append(new_metabolite)
                    total_metabolites_produced += 1
                    # Mark glucose as consumed for this step
                    for index in consumed_indices:
                        glucose_list[index] = None
                
                # Handle division and death
                action, result = gillespie_step(bacteria, glucose_list, glucose_grid, metabolite_list, metabolite_grid, space_size)
                
                if action == 'divide':
                    if result: new_bact.append(result)
                    survivors.append(bacteria)
                elif action == 'die_by_toxicity':
                    toxicity_deaths += 1
                    # Don't append to survivors
                elif action == 'die':
                    # Don't append to survivors
                    pass
                else: # 'idle'
                    survivors.append(bacteria)

        ecoli_list = survivors + new_bact
        
        # Remove consumed glucose from the list
        glucose_list = [g for g in glucose_list if g is not None]
        
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
            glucose_list = add_glucose_from_top(glucose_list, adjusted_feed_amount, config)
        
        # Save state for visualization
        states.append((copy.deepcopy(ecoli_list), 
                      copy.deepcopy(glucose_list), 
                      copy.deepcopy(metabolite_list),
                      total_metabolites_produced,
                      drain_system.collected_metabolites if drain_system else 0,
                      toxicity_deaths))
        
        # Check for simulation end condition
        if not ecoli_list:
            print(f"\n\nSimulation stopped at step {step + 1}: no bacteria remaining.")
            break

    return states, space_size, drain_system, fluid

# ====================
# --- Data Export ---
# ====================

def _calculate_population_stats(ecolis, glucoses, metabolites, step, prev_bacteria_count):
    """Helper to calculate population-based statistics."""
    n_bacteria = len(ecolis)
    growth_rate = (n_bacteria / prev_bacteria_count) if prev_bacteria_count > 0 else 1.0
    return {
        'Bacteria_Count': n_bacteria,
        'Glucose_Count': len(glucoses),
        'Metabolite_Count': len(metabolites),
        'Growth_Rate': f"{growth_rate:.3f}",
    }

def _calculate_atp_stats(ecolis):
    """Helper to calculate ATP-related statistics."""
    if not ecolis:
        return {
            'Average_ATP': "0.00", 'Min_ATP': "0.00", 'Max_ATP': "0.00",
            'Total_ATP': "0.00", 'ATP_StdDev': "0.00"
        }
    atp_values = [e.ATP for e in ecolis]
    return {
        'Average_ATP': f"{np.mean(atp_values):.2f}",
        'Min_ATP': f"{min(atp_values):.2f}",
        'Max_ATP': f"{max(atp_values):.2f}",
        'Total_ATP': f"{sum(atp_values):.2f}",
        'ATP_StdDev': f"{np.std(atp_values):.2f}",
    }

def _calculate_metabolite_stats(ecolis, metabolites, total_metabolites_produced, collected_count):
    """Helper to calculate metabolite-related statistics."""
    metabolite_production = [e.produced_metabolite for e in ecolis]
    avg_metabolites_per_bacteria = np.mean(metabolite_production) if ecolis else 0
    metabolite_concentration = [m.concentration for m in metabolites] if metabolites else [0]
    
    return {
        'Total_Metabolites_Produced': total_metabolites_produced,
        'Collected_Metabolites': collected_count,
        'Avg_Metabolites_Per_Bacteria': f"{avg_metabolites_per_bacteria:.2f}",
        'Avg_Metabolite_Concentration': f"{np.mean(metabolite_concentration):.3f}",
    }

def _calculate_spatial_stats(particles):
    """Helper to calculate spatial statistics for a list of particles."""
    if not particles:
        return {'Avg_X': "0.000000", 'Avg_Y': "0.000000", 'Spatial_Std': "0.000000"}
        
    x_pos = [p.x for p in particles]
    y_pos = [p.y for p in particles]
    return {
        'Avg_X': f"{np.mean(x_pos):.6f}",
        'Avg_Y': f"{np.mean(y_pos):.6f}",
        'Spatial_Std': f"{np.std(x_pos) + np.std(y_pos):.6f}",
    }

def _calculate_efficiency_stats(n_bacteria, n_glucose, total_atp, total_metabolites_produced):
    """Helper to calculate resource efficiency metrics."""
    if n_bacteria == 0:
        return {
            'Glucose_Per_Bacteria': "0.00", 'ATP_Per_Glucose': "0.00",
            'Metabolites_Per_Glucose': "0.00", 'ATP_Per_Bacteria': "0.00",
            'Metabolites_Per_Bacteria': "0.00"
        }
        
    return {
        'Glucose_Per_Bacteria': f"{n_glucose / n_bacteria:.2f}",
        'ATP_Per_Glucose': f"{total_atp / n_glucose:.2f}" if n_glucose > 0 else "0.00",
        'Metabolites_Per_Glucose': f"{total_metabolites_produced / n_glucose:.2f}" if n_glucose > 0 else "0.00",
        'ATP_Per_Bacteria': f"{total_atp / n_bacteria:.2f}",
        'Metabolites_Per_Bacteria': f"{total_metabolites_produced / n_bacteria:.2f}"
    }

def export_detailed_data(states, filename="simulation_data.csv"):
    """
    Export detailed simulation data to CSV with comprehensive statistics.
    This function is now a high-level coordinator using helper functions.
    """
    
    statistics = []
    for step, (ecolis, glucoses, metabolites, total_metabolites_produced, collected_count, toxicity_deaths) in enumerate(states):
        if not ecolis and step > 0: # Stop if population has died out
            continue
            
        # --- Calculate all statistics using helpers ---
        prev_count = len(states[step-1][0]) if step > 0 else len(ecolis)
        pop_stats = _calculate_population_stats(ecolis, glucoses, metabolites, step, prev_count)
        atp_stats = _calculate_atp_stats(ecolis)
        metabolite_stats = _calculate_metabolite_stats(ecolis, metabolites, total_metabolites_produced, collected_count)
        
        bacteria_spatial = _calculate_spatial_stats(ecolis)
        glucose_spatial = _calculate_spatial_stats(glucoses)
        
        total_atp = float(atp_stats['Total_ATP'])
        efficiency_stats = _calculate_efficiency_stats(pop_stats['Bacteria_Count'], pop_stats['Glucose_Count'], total_atp, total_metabolites_produced)
        
        # --- Combine all statistics into a single record ---
        stats = {
            'Step': step,
            'Time_stamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f"),
            **pop_stats,
            **atp_stats,
            **metabolite_stats,
            'Toxicity_Deaths': toxicity_deaths,
            'Avg_Bacteria_X': bacteria_spatial['Avg_X'],
            'Avg_Bacteria_Y': bacteria_spatial['Avg_Y'],
            'Bacteria_Spatial_Std': bacteria_spatial['Spatial_Std'],
            'Avg_Glucose_X': glucose_spatial['Avg_X'],
            'Avg_Glucose_Y': glucose_spatial['Avg_Y'],
            'Glucose_Spatial_Std': glucose_spatial['Spatial_Std'],
            **efficiency_stats
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
            'Total_Metabolites_Produced', 'Collected_Metabolites', 'Avg_Metabolites_Per_Bacteria', 'Avg_Metabolite_Concentration',
            # New Toxicity Metric
            'Toxicity_Deaths',
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
    total_metabolites = [s[3] for s in states]
    avg_atp = [np.mean([e.ATP for e in s[0]]) if s[0] else 0 for s in states]

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

def launch_dashboard(states, space_size, fluid_field, drain_system, scenario_num, scenario_name):
    """Initializes and runs the Dash dashboard in a self-contained scope."""

    class SimulationDashboard:
        def __init__(self, states, space_size, fluid_field, drain_system, scenario_num, scenario_name):
            self.states = states
            self.space_size = space_size
            self.fluid_field = fluid_field
            self.drain_system = drain_system
            self.scenario_num = scenario_num
            self.scenario_name = scenario_name
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
                    'toxic_concentration': list(self.toxic_concentration_history),
                    'total_acetate': list(self.total_acetate_history),
                    'collected_acetate': list(self.collected_acetate_history),
                    'toxicity_deaths': list(self.toxicity_deaths_history)
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
            self.total_acetate_history = deque(maxlen=self.history_length)
            self.collected_acetate_history = deque(maxlen=self.history_length)
            self.toxicity_deaths_history = deque(maxlen=self.history_length)

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
                        html.H1("E. coli Bioreactor Simulation",
                               className="text-center mt-4 mb-2 text-primary"),
                        html.P(f"Scenario {self.scenario_num}: {self.scenario_name.replace('_', ' ').title()}",
                               className="text-center text-muted lead mb-4")
                    ])
                ]),
                
                # Main content
                dbc.Row([
                    # Left column: Simulation view and controls
                    dbc.Col([
                        dbc.Card([
                            dbc.CardHeader("Live Simulation Environment"),
                            dbc.CardBody([
                                dcc.Graph(id='simulation-graph', style={'height': '60vh'})
                            ])
                        ], className="mb-4 shadow-sm"),
                        
                        dbc.Card([
                            dbc.CardHeader("Simulation Controls"),
                            dbc.CardBody([
                                dbc.Row([
                                    dbc.Col([
                                        dbc.ButtonGroup([
                                            dbc.Button("▶ Play", id="play-button", color="success", className="me-2"),
                                            dbc.Button("⏸ Pause", id="pause-button", color="warning", className="me-2"),
                                            dbc.Button("↺ Reset", id="reset-button", color="danger")
                                        ], className="mb-3")
                                    ], width='auto'),
                                    dbc.Col([
                                        html.Div(id='current-step-display', className="text-muted mt-2")
                                    ])
                                ]),
                                html.Label("Animation Speed", className="fw-bold"),
                                dcc.Slider(
                                    id='speed-slider',
                                    min=0.1, max=2.0, step=0.1, value=1.0,
                                    marks={i/10: str(i/10) for i in range(2, 21, 2)},
                                    className="mb-2"
                                )
                            ])
                        ], className="mb-4 shadow-sm"),
                        dbc.Card([
                            dbc.CardHeader("Live Statistics"),
                            dbc.CardBody([
                                html.Div(id='statistics-display', className="p-3")
                            ])
                        ], className="shadow-sm")
                    ], width=7),
                    
                    # Right column: Time series and statistics in a single column
                    dbc.Col([
                        dbc.Card([
                            dbc.CardHeader("Time Series KPIs"),
                            dbc.CardBody([
                                dcc.Graph(id='time-series-graph')
                            ])
                        ], className="mb-4 shadow-sm")
                    ], width=5)
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
                ecoli_list, glucose_list, metabolite_list, total_metabolites_produced, collected_metabolites, toxicity_deaths = self.states[self.current_step]
                
                # Update histories
                self.time_history.append(self.current_step)
                self.metabolites_history.append(len(metabolite_list))
                self.population_history.append(len(ecoli_list))
                self.atp_history.append(np.mean([e.ATP for e in ecoli_list]) if ecoli_list else 0)
                self.excreted_metabolites_history.append(total_metabolites_produced)
                self.toxic_concentration_history.append(self.calculate_average_toxicity(ecoli_list, metabolite_list))
                self.total_acetate_history.append(sum(m.concentration for m in metabolite_list))
                self.collected_acetate_history.append(collected_metabolites)
                self.toxicity_deaths_history.append(toxicity_deaths)
                
                # Create views
                sim_fig = self.create_simulation_view(ecoli_list, glucose_list, metabolite_list)
                ts_fig = self.create_time_series_view()
                stats = self.create_statistics_display(ecoli_list, glucose_list, metabolite_list, total_metabolites_produced, collected_metabolites, toxicity_deaths)
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
                        'toxic_concentration': list(self.toxic_concentration_history),
                        'total_acetate': list(self.total_acetate_history),
                        'collected_acetate': list(self.collected_acetate_history),
                        'toxicity_deaths': list(self.toxicity_deaths_history)
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
                    name='Acetate',
                    marker=dict(
                        size=4,
                        color=[m.concentration for m in metabolite_list],
                        colorscale='RdPu',  # Red-Purple scale for toxicity
                        opacity=0.4,
                        showscale=False
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
                    color=[e.ATP for e in ecoli_list],
                    colorscale='Viridis',
                    opacity=0.9,
                    cmin=0,
                    cmax=(ecoli_list[0].config.ECOLI_DIVISION_THRESHOLD_ATP * 1.1) if ecoli_list else 220,
                    line=dict(width=1, color='rgba(0, 0, 0, 0.7)'), # Add border to bacteria
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
                template="plotly_dark", # Use a lighter theme for better contrast
                margin=dict(l=40, r=140, t=40, b=40),  # Increased right margin for colorbar
                xaxis=dict(
                    title='X Position (micrometers)',
                    showgrid=True,
                    gridcolor='rgba(255, 255, 255, 0.1)',
                    titlefont=dict(size=14),  # Larger axis title
                    tickfont=dict(size=12)  # Larger tick labels
                ),
                yaxis=dict(
                    title='Y Position (micrometers)',
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
                rows=8, cols=1,
                subplot_titles=(
                    'Active Acetate in Environment',
                    'Bacterial Population Dynamics',
                    'Average ATP Level per Bacterium',
                    'Cumulative Acetate Production',
                    'Avg. Local Acetate Conc. Experienced by Bacteria',
                    'Total Acetate Concentration in Reactor',
                    'Collected Acetate Over Time',
                    'Deaths Caused by Acetate Toxicity'
                ),
                vertical_spacing=0.08
            )
            
            # Add time series traces with custom styling and clear names
            fig.add_trace(
                go.Scatter(
                    x=list(self.time_history), 
                    y=list(self.metabolites_history),
                    name='Active Acetate Count',
                    line=dict(color='#28a745', width=2),
                    showlegend=True
                ),
                row=1, col=1
            )
            
            fig.add_trace(
                go.Scatter(
                    x=list(self.time_history), 
                    y=list(self.population_history),
                    name='E. coli Population',
                    line=dict(color='#dc3545', width=2),
                    showlegend=True
                ),
                row=2, col=1
            )
            
            fig.add_trace(
                go.Scatter(
                    x=list(self.time_history), 
                    y=list(self.atp_history),
                    name='Average ATP Level',
                    line=dict(color='#ffc107', width=2),
                    showlegend=True
                ),
                row=3, col=1
            )
            
            fig.add_trace(
                go.Scatter(
                    x=list(self.time_history), 
                    y=list(self.excreted_metabolites_history),
                    name='Total Acetate Produced',
                    line=dict(color='#17a2b8', width=2),
                    showlegend=True
                ),
                row=4, col=1
            )
            
            fig.add_trace(
                go.Scatter(
                    x=list(self.time_history), 
                    y=list(self.toxic_concentration_history),
                    name='Average Toxic Concentration',
                    line=dict(color='#e83e8c', width=2),
                    showlegend=True
                ),
                row=5, col=1
            )
            
            fig.add_trace(
                go.Scatter(
                    x=list(self.time_history),
                    y=list(self.total_acetate_history),
                    name='Total Acetate Concentration',
                    line=dict(color='#6f42c1', width=2),
                    showlegend=True
                ),
                row=6, col=1
            )
            
            fig.add_trace(
                go.Scatter(
                    x=list(self.time_history),
                    y=list(self.collected_acetate_history),
                    name='Collected Acetate',
                    line=dict(color='#fd7e14', width=2),
                    showlegend=True
                ),
                row=7, col=1
            )

            fig.add_trace(
                go.Scatter(
                    x=list(self.time_history),
                    y=list(self.toxicity_deaths_history),
                    name='Toxicity Deaths',
                    line=dict(color='#6610f2', width=2), # A new color
                    showlegend=True
                ),
                row=8, col=1
            )
            
            # Update layout with improved styling
            fig.update_layout(
                height=900,  # Adjusted height for new layout
                showlegend=True, # Legend is now redundant with titles
                template="plotly_dark",
                margin=dict(l=60, r=40, t=80, b=50),
                legend=dict(
                    orientation="h",
                    yanchor="bottom",
                    y=-0.15, # Position below plots
                    xanchor="center",
                    x=0.5
                )
            )
            
            # Update axes with improved labels
            for i in range(1, 9):
                fig.update_xaxes(
                    title_text="Time Steps",
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='rgba(255, 255, 255, 0.1)',
                    row=i, col=1
                )
            
            # Add specific y-axis labels
            y_labels = [
                "Acetate Particle Count",
                "Bacterial Cell Count",
                "Average ATP Units",
                "Total Acetate Molecules",
                "Avg. Local Concentration",
                "Total Concentration (mM)",
                "Collected Acetate Count",
                "Toxicity Death Count"
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

        def create_statistics_display(self, ecoli_list, glucose_list, metabolite_list, total_metabolites_produced, collected_metabolites, toxicity_deaths):
            """Create the statistics display panel"""
            # Calculate average toxic concentration
            avg_toxic_conc = self.calculate_average_toxicity(ecoli_list, metabolite_list)
            avg_atp = np.mean([e.ATP for e in ecoli_list]) if ecoli_list else 0
            
            def create_stat_card(title, value, color):
                return dbc.Col([
                    html.H4(value, className=f"text-{color}"),
                    html.P(title, className="text-muted")
                ], className="text-center mb-3")

            stats = dbc.Row([
                create_stat_card("Live Bacteria", len(ecoli_list), "success"),
                create_stat_card("Available Glucose", len(glucose_list), "info"),
                create_stat_card("Active Acetate", len(metabolite_list), "primary"),
                create_stat_card("Average ATP Level", f"{avg_atp:.2f}", "warning"),
                create_stat_card("Total Acetate Produced", f"{total_metabolites_produced}", "secondary"),
                create_stat_card("Collected Acetate", f"{collected_metabolites}", "dark"),
                create_stat_card("Avg. Local Acetate", f"{avg_toxic_conc:.2f}", "danger"),
                create_stat_card("Toxicity Deaths", f"{toxicity_deaths}", "danger")
            ], className="g-3")
            
            return stats
    
    # Initialize and run dashboard
    dashboard = SimulationDashboard(states, space_size, fluid_field=fluid_field, drain_system=drain_system, scenario_num=choice, scenario_name=selected_scenario['name'])
    dashboard.run(debug=False, port=8050)

# ====================
# --- Main ---
# ====================

def print_welcome_banner():
    """Prints a friendly welcome banner for the user."""
    banner = """
    ================================================================
    |                                                              |
    |          Welcome to the E. coli Bioreactor Simulator         |
    |                                                              |
    |      A dynamic simulation of bacterial growth, metabolism,   |
    |                   and environmental interaction.             |
    |                                                              |
    ================================================================
    """
    print(banner)

if __name__ == '__main__':
    # Print the welcome banner
    print_welcome_banner()

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
            'name': 'no_drain',
            'description': 'No particle recovery or drainage system',
            'n_ecoli': 50,
            'n_glucose': 1000,
            'max_steps': 500,
            'use_drain': False
        },
        '3': {
            'name': 'glucose_rich',
            'description': 'High initial and frequent glucose feed',
            'n_ecoli': 50,
            'n_glucose': 3000,              
            'glucose_feed_interval': 10,  
            'glucose_feed_amount': 200,   
            'max_steps': 500
        },
        '4': {
            'name': 'glucose_poor',
            'description': 'Low glucose and rare feeding',
            'n_ecoli': 50,
            'n_glucose': 200,              
            'glucose_feed_interval': 50,   
            'glucose_feed_amount': 20, 
            'max_steps': 500
        }
    }
    
    # Print available scenarios in a user-friendly format
    print("\nPlease select a simulation scenario from the options below:\n")
    print("-" * 60)
    for key, scenario in scenarios.items():
        print(f"  [{key}] {scenario['name'].replace('_', ' ').title()}")
        print(f"      Description: {scenario['description']}")
        print(f"      Key Parameters:")
        for p_key, p_val in scenario.items():
            if p_key not in ['name', 'description']:
                print(f"        - {p_key.replace('_', ' ').title()}: {p_val}")
        print("-" * 60)
    
    # Get user choice
    while True:
        choice = input("\nEnter the number of the scenario you want to run (1-4): ").strip()
        if choice in scenarios:
            break
        print("\nInvalid choice. Please enter a number from 1 to 4.")
    
    # Get selected scenario
    selected_scenario = scenarios[choice]
    print(f"\n✅ You have selected '{selected_scenario['name'].replace('_', ' ').title()}'.")
    print("      Initializing simulation with the following parameters:")
    for k, v in selected_scenario.items():
        if k not in ['name', 'description']:
            print(f"        - {k.replace('_', ' ').title()}: {v}")
    
    # Determine the output directory in a portable way
    try:
        # Ideal case: Get the directory of the current script
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        # Fallback for environments where __file__ is not defined (e.g., notebooks, REPL)
        script_dir = os.getcwd()
        print(f"\nWarning: Could not determine script path. Saving files to the current working directory: {script_dir}")
    
    # Create a configuration instance
    config = SimulationConfig()
    
    # Get space size from scenario or use default, and update config
    space_size = selected_scenario.get('space_size', config.SPACE_SIZE)
    config.SPACE_SIZE = space_size
    
    # Run simulation - The fluid object is now returned from here
    states, space_size, drain_system, fluid = simulate(
        n_ecoli=selected_scenario['n_ecoli'],
        n_glucose=selected_scenario['n_glucose'],
        space_size=space_size,
        max_steps=selected_scenario.get('max_steps', 250),
        glucose_feed_interval=selected_scenario.get('glucose_feed_interval', 20),
        glucose_feed_amount=selected_scenario.get('glucose_feed_amount', 50),
        use_drain=selected_scenario.get('use_drain', True),
        config=config
    )
    
    # Export simulation data to CSV
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_filename = os.path.join(script_dir, f"simulation_data_{timestamp}.csv")
    export_detailed_data(states, filename=csv_filename)
    print(f"\n📊 Detailed simulation data has been exported to: {csv_filename}")

    #export simulation data to PNG
    png_filename = os.path.join(script_dir, "simulation_summary.png")
    export_summary_image(states, filename=png_filename)
    print(f"📈 A summary image of the simulation has been saved as: {png_filename}")
    
    # Initialize and run dashboard
    if not drain_system:
        # If the drain was disabled, we still need a dummy object for the dashboard
        # This object won't do anything, but prevents dashboard components from crashing.
        drain_system = DrainSystem(space_size=space_size, config=config)

    print("\n🚀 Launching interactive dashboard...")
    print("   Please open your web browser and navigate to the link printed below.")
    print("   (It may take a few seconds for the server to start)")

    launch_dashboard(states, space_size, fluid, drain_system, scenario_num=choice, scenario_name=selected_scenario['name'])
