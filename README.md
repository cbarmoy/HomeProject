# 🦠 Bioreactor Simulation

Welcome to the Bioreactor Simulation! This interactive tool, built with Python and Dash, models the behavior of *E. coli* in a dynamic microenvironment. Watch as bacteria consume glucose, produce acetate, and respond to environmental factors like fluid dynamics and metabolite toxicity.

## 🚀 Prerequisites

Before we dive into the microscopic mosh pit, you'll need:

1. **Python** (Version 3.8 or higher)
   - Download from [Python's official website](https://www.python.org/downloads/)
   - Don't forget to check "Add Python to PATH" (or your computer will give you the cold shoulder 🥶)

2. **Spyder** (Your friendly neighborhood IDE)
   - Get the whole squad through Anaconda: [Download here](https://www.anaconda.com/products/distribution)
   - It's like Discord for code, but with more graphs and fewer emojis 📊

## 🎮 Setup Instructions

When you run the script, you will be prompted to choose from several predefined scenarios. After selecting a scenario, the simulation will generate the necessary data and launch a web-based dashboard for visualization.

### Method 1: The Easy Way (Spyder) 
#### AKA "I Just Want to Watch Bacteria Dance" Edition


1. **Open the Anaconda Prompt** (on Windows) or a terminal (on macOS/Linux).
2. Run the following command to install the required packages:

     ```
     pip install numpy dash plotly pandas dash-bootstrap-components
     ```

3. Open `Bioreactor_Simulation.py` in Spyder.
4. Hit that green play button! 🎵
5. Follow the prompts in the console to select a scenario.
6. The dashboard will automatically open in your web browser.

### Method 2: The Terminal Way 
#### AKA "I'm a Command Line Warrior" Edition

1. Open your terminal (Command Prompt for Windows warriors, Terminal for Mac/Linux legends)
2. Navigate to your project directory:
   ```bash
   cd path/to/your/project/folder
   ```
3. Create a virtual environment (recommended):
   ```bash
   python -m venv .venv
   ```
4. Activate the virtual environment:
   - Windows:
     ```bash
     .venv\Scripts\activate
     ```
   - macOS/Linux:
     ```bash
     source .venv/bin/activate
     ```
5. Install the required packages:
   ```bash
   pip install -r requirements.txt
   ```
6. Start the simulation:
   ```bash
   python Bioreactor_Simulation.py
   ```
7. Follow the prompts in the terminal to select a scenario. The dashboard will then launch in your browser.

## 🆘 Troubleshooting (When the Party Goes Wrong)

1. **"Module not found" error**
   - It looks like a required package is missing.
   - Run `pip install -r requirements.txt` to install all necessary dependencies.

2. **"Port already in use" error**
   - Another application is using port 8050.
   - You can specify a different port when launching the dashboard (this requires a minor code modification).

3. **Dashboard not showing**
   - The simulation first runs in the terminal to generate data. Once complete, it will launch the dashboard.
   - Check the console output for the URL (usually http://127.0.0.1:8050) and open it in your browser.

4. **Slow performance**
   - The simulation can be computationally intensive.
   - Consider choosing a scenario with fewer particles or a smaller number of steps.

## 🎨 Features

- **Multi-particle System**:
  - *E. coli* bacteria with ATP-based metabolism and chemotaxis.
  - Glucose particles serving as the primary substrate.
  - Acetate (metabolite) particles that diffuse and can become toxic.
- **Fluid Dynamics**:
  - A simulated flow field affects particle movement, mimicking a bioreactor environment.
  - A drain system removes particles from the bottom and can recycle them.
- **Advanced Biological Models**:
  - ATP-based energy management for bacterial growth, division, and survival.
  - Michaelis-Menten-like kinetics for glucose consumption.
  - Metabolite production (acetate) and its toxic effects on the bacterial population.
- **Interactive Dashboard**:
  - Real-time visualization of particle positions.
  - Dynamic plots for population counts, ATP levels, and metabolite concentrations.
  - Controls to play, pause, and adjust the simulation speed.

## 🎯 Available Scenarios

When you run the script, you can choose from one of the following scenarios:

1. **Baseline** – The Standard Run 🎈  
   - A balanced simulation with standard parameters.  
   - `n_ecoli: 50`, `n_glucose: 1000`, `max_steps: 500`

2. **No Drain** – The Closed System 🚫🧹  
   - All particles remain in the simulation, leading to acetate accumulation.  
   - `use_drain: False`

3. **Glucose Rich** – The High-Yield Experiment 🍰  
   - An abundance of glucose to promote rapid growth.  
   - `n_glucose: 3000`, `glucose_feed_interval: 10`, `glucose_feed_amount: 200`

4. **Glucose Poor** – The Survival Challenge 🍃  
   - Limited glucose availability, testing the bacteria's resilience.  
   - `n_glucose: 200`, `glucose_feed_interval: 50`, `glucose_feed_amount: 20`

## 🎮 Running the Show

1. Pick your party scenario (1-4)
2. Watch the magic happen in your browser
3. Use the dashboard to:
   - Play/Pause (for dramatic effect ⏯️)
   - Reset (when things get too wild 🔄)
   - Control speed (from chill vibes to chaos 🌪️)

## 📊 Understanding What You're Seeing

### Dashboard Elements
- **Simulation View**:
    - Green dots: *E. coli* bacteria
    - Blue dots: Glucose particles
    - Purple dots: Acetate (metabolite) particles
    - Gray areas: The drain system at the bottom.
- **Time-Series Plots**:
    - Track key metrics over time, including population counts, ATP levels, and acetate concentration.

### CSV Data Output
The simulation automatically saves detailed data from each run:
- File name: `simulation_data_YYYYMMDD_HHMMSS.csv`
- This file contains a step-by-step log of key simulation parameters, perfect for further analysis.

## 🧬 Code Breakdown (For the Curious Minds)

### 🏗️ Main Components

The code is structured around a central configuration object and several key classes that model the different components of the bioreactor.

#### `SimulationConfig` Class
A `dataclass` that holds all the simulation parameters. This makes it easy to manage and adjust settings from a single location.

#### `StreamField` Class 🌊
- Manages the fluid dynamics of the environment.
- Simulates a gentle downward flow with upward currents and drain effects.

#### `SpatialGrid` Class 📍
- Implements spatial hashing to speed up neighbor detection.
- Crucial for efficiently calculating interactions between particles.

#### Particle Classes 🎯

1. **`Ecoli` Class** 🦠
- Models the bacteria. They can:
  - Move via chemotaxis towards glucose.
  - Consume glucose to gain ATP.
  - Produce acetate as a metabolite.
  - Divide if they have sufficient ATP.
  - Die from low ATP or high acetate toxicity.

2. **`Glucose` Class** 🍪
- The substrate particles.
- They move based on diffusion and fluid dynamics.

3. **`Metabolite` Class** 🧪
- Represents acetate, the waste product.
- Diffuses through the environment and can negatively impact bacteria at high concentrations.

#### `DrainSystem` Class 🚰
- Manages the collection and recycling of particles that reach the bottom of the reactor.

### 🎮 Simulation and Visualization

#### `simulate()` Function
- The core simulation engine that runs the step-by-step calculations before the dashboard is launched.

#### `SimulationDashboard` Class 🖥️
- An interactive dashboard built with Dash and Plotly.
- Provides visualization and controls for exploring the simulation results.

### 🔬 The Science Behind It

1. **Movement and Physics** 🏃‍♂️
   - **Chemotaxis**: Bacteria move towards higher concentrations of glucose.
   - **Fluid Dynamics**: Particle movement is influenced by a simulated flow field.
   - **Diffusion**: Glucose and acetate particles exhibit random, diffusion-like movement.

2. **Metabolism and Growth** ⚡
   - **ATP-Based Energy**: Bacterial actions like movement and division are fueled by ATP, which is gained from consuming glucose.
   - **Toxicity Model**: High concentrations of acetate create a toxic environment that drains ATP and can lead to cell death.

### 📁 File Structure

```
HomeProject/
│
├── Bioreactor_Simulation.py # The main simulation and dashboard script. 🎉
├── README.md                # You are here! 📍
└── requirements.txt         # A list of the required Python packages. 🛒
```

### 🎛️ Key Parameters You Can Tweak

All key parameters are consolidated in the `SimulationConfig` class at the top of `Bioreactor_Simulation.py`.

```python
# In Bioreactor_Simulation.py

@dataclass
class SimulationConfig:
    # E. coli Parameters
    ECOLI_INITIAL_ATP: float = 200.0
    ECOLI_GLUCOSE_SENSITIVITY: float = 0.5
    ECOLI_BASE_GROWTH_RATE: float = 0.12

    # Acetate (Metabolite) Parameters
    ACETATE_TOXICITY_THRESHOLD: float = 8.0
    ACETATE_LETHAL_CONCENTRATION: float = 60.0

    # ... and many more!
```

Feel free to experiment with these values to see how they affect the simulation outcome.

## 📜 License

This project is licensed under the MIT License - because even bacteria need legal protection! 🤝

Remember: No bacteria were harmed in the making of this simulation. They're just pixels having a good time! 🦠✨

## 📞 Need Help?

Did your bacteria throw an unexpected party? Code acting sus? Don't panic! 

### Contact the Developer
- **GitHub**: [cbarmoy](https://github.com/cbarmoy)
- **Email**: come. barmoy1@supbiotech.fr

Feel free to reach out if:
- Your bacteria are misbehaving 🦠
- You found a bug (the software kind, not the bacterial kind) 🐛
- You have ideas for new features 💡
- You just want to share your bacterial party stories! 🎉

I'm always happy to help fellow bacterial party organizers! 🚀 
