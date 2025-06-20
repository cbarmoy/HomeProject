# 🦠 Bacterial Growth Simulation

Welcome to the most epic bacteria party simulator you'll ever see! Watch as your tiny bacterial friends navigate through life, eat glucose, and try not to get too tipsy on their own metabolites. 🎉

## 🚀 Prerequisites

Before we dive into the microscopic mosh pit, you'll need:

1. **Python** (Version 3.8 or higher)
   - Download from [Python's official website](https://www.python.org/downloads/)
   - Don't forget to check "Add Python to PATH" (or your computer will give you the cold shoulder 🥶)

2. **Spyder** (Your friendly neighborhood IDE)
   - Get the whole squad through Anaconda: [Download here](https://www.anaconda.com/products/distribution)
   - It's like Discord for code, but with more graphs and fewer emojis 📊

## 🎮 Setup Instructions

### Method 1: The Easy Way (Spyder) 
#### AKA "I Just Want to Watch Bacteria Dance" Edition


1. **Open the Anaconda Prompt** (on Windows) or a terminal (on macOS/Linux).
2. Run one of the following commands:

     ```
     pip install dash
     pip install plotly
     pip install dash-bootstrap-components
     ```

3. Open "Final Home Proj.py" in spyder (your bacteria's party playlist)
4. Hit that green play button like it's a Discord call! 🎵

### Method 2: The Terminal Way 
#### AKA "I'm a Command Line Warrior" Edition

1. Open your terminal (Command Prompt for Windows warriors, Terminal for Mac/Linux legends)
2. Navigate to your bacteria's dance floor:
   ```bash
   cd path/to/project/folder
   ```
3. Create a cozy virtual environment:
   ```bash
   python -m venv .venv
   ```
4. Activate your bacteria's VIP room:
   - Windows party:
     ```bash
     .venv\Scripts\activate
     ```
   - Mac/Linux rave:
     ```bash
     source .venv/bin/activate
     ```
5. Install the party supplies:
   ```bash
   pip install numpy dash plotly pandas dash-bootstrap-components
   ```
6. Start the party:
   ```bash
   python "Final Home Proj.py"
   ```

## 🆘 Troubleshooting (When the Party Goes Wrong)

1. **"Module not found" error**
   - Looks like someone forgot to bring snacks to the party
   - Try `pip install [package_name]` to feed the hungry Python 🐍

2. **"Port already in use" error**
   - Another party is already happening on port 8050
   - Time to crash that other party or find a new venue!

3. **Dashboard not showing**
   - Your bacteria are shy! Visit them at: http://127.0.0.1:8050
   - Make sure no other app is hogging port 8050 (looking at you, Spotify 👀)

4. **Slow performance**
   - Too many bacteria at the party!
   - Try reducing the crowd size (your CPU will thank you)

## 🎨 Features

- **Multi-particle System**:
  - E. coli bacteria with ATP-based metabolism
  - Glucose particles with momentum-based movement
  - Metabolites with diffusion and concentration decay
  - Fluid dynamics affecting particle movement

- **Advanced Physics**:
  - Realistic particle movement with momentum
  - Fluid dynamics with streamlines
  - Gravity effects and boundary interactions
  - Brownian motion and diffusion

- **Biological Mechanisms**:
  - ATP-based bacterial metabolism
  - Glucose consumption and metabolite production
  - Toxicity effects from metabolite concentration
  - Bacterial growth and division

- **Environmental Features**:
  - Drain system for particle collection and recycling
  - Dynamic glucose feeding system
  - Spatial partitioning for efficient neighbor detection
  - Customizable environment parameters

## 🎯 Available Scenarios

1. **Baseline** – The Casual Friday Party 🎈  
   - Just your standard bacterial get-together  
   - Parameters: `n_ecoli: 50`, `n_glucose: 1000`, `max_steps: 500`

2. **No Drain** – The Hoarder’s Paradise 🚫🧹  
   - No cleanup crew in sight — what goes in stays in  
   - Parameters: `n_ecoli: 50`, `n_glucose: 1000`, `max_steps: 500`, `use_drain: False`

3. **Glucose Rich** – The All-You-Can-Eat Buffet 🍰  
   - Endless snacks and frequent refills — pure bacterial bliss  
   - Parameters: `n_ecoli: 50`, `n_glucose: 3000`, `glucose_feed_interval: 10`, `glucose_feed_amount: 200`, `max_steps: 500`

4. **Glucose Poor** – The Fasting Retreat 🍃  
   - Minimal snacks, long wait times — survival of the chillest  
   - Parameters: `n_ecoli: 50`, `n_glucose: 200`, `glucose_feed_interval: 50`, `glucose_feed_amount: 20`, `max_steps: 500`

## 🎮 Running the Show

1. Pick your party scenario (1-4)
2. Watch the magic happen in your browser
3. Use the dashboard to:
   - Play/Pause (for dramatic effect ⏯️)
   - Reset (when things get too wild 🔄)
   - Control speed (from chill vibes to chaos 🌪️)

## 📊 Understanding What You're Seeing

### Dashboard Elements
- Blue dots: Glucose (the snacks)
- Green dots: Your bacterial party people
- Purple dots: Metabolites (what happens after too many snacks)
- Gray bars: The bouncers (drain system)

### CSV Data Output
All party statistics are saved for posterity:
- File name: `simulation_data_YYYYMMDD_HHMMSS.csv`
- Perfect for spreadsheet enthusiasts and data nerds 🤓

## 🧬 Code Breakdown (For the Curious Minds)

### 🏗️ Main Components

#### StreamField Class 🌊
```python
class StreamField:
    """Your bacteria's swimming pool!"""
```
- Creates the fluid environment where everything happens
- Manages flow patterns and drain positions
- Think of it as the bacteria's water park with special currents

#### SpatialGrid Class 📍
```python
class SpatialGrid:
    """The bacteria's GPS system"""
```
- Helps bacteria find their friends and food efficiently
- Divides space into smaller chunks (like a Discord server with multiple channels)
- Makes sure your CPU doesn't have a meltdown calculating distances

#### Particle Classes 🎯

1. **Ecoli Class** 🦠
```python
class Ecoli:
    """The main party people!"""
```
- Your bacterial protagonists
- They can:
  - Move around looking for food (chemotaxis)
  - Eat glucose for ATP (energy drinks!)
  - Produce metabolites (party leftovers)
  - Divide when they're happy and well-fed
  - Die if they're too hungry or surrounded by too many metabolites

2. **Glucose Class** 🍪
```python
class Glucose:
    """The snacks"""
```
- The food particles
- Features:
  - Floats around with realistic physics
  - Gets eaten by bacteria
  - Respawns at the top (like Discord messages, but tastier)

3. **Metabolite Class** 🧪
```python
class Metabolite:
    """The party aftermath"""
```
- Waste products from bacterial metabolism
- Can become toxic if too concentrated
- Gets collected by the drain system

#### DrainSystem Class 🚰
```python
class DrainSystem:
    """The cleanup crew"""
```
- Collects and recycles particles
- Manages the party cleanup
- Keeps track of what's been collected

### 🎮 Simulation Control

#### Main Simulation Loop
```python
def simulate():
    """The party manager"""
```
- Controls the whole simulation
- Updates positions and states
- Handles:
  - Bacterial movement and feeding
  - Glucose distribution
  - Metabolite production
  - Particle recycling
  - Statistics tracking

#### Dashboard (The UI) 🖥️
```python
class SimulationDashboard:
    """Your control room"""
```
- Built with Dash and Plotly
- Shows:
  - Real-time particle visualization
  - Population graphs
  - ATP levels
  - Metabolite concentrations
  - Interactive controls

### 🔬 The Science Behind It

1. **Movement Physics** 🏃‍♂️
   - Uses momentum and fluid dynamics
   - Includes Brownian motion (random wiggling)
   - Simulates chemotaxis (bacteria following food)

2. **Metabolism System** ⚡
   - ATP-based energy management
   - Glucose consumption mechanics
   - Metabolite production and toxicity

3. **Population Dynamics** 📈
   - Growth based on available resources
   - Death from starvation or toxicity
   - Spatial distribution patterns

### 📁 File Structure

```
HomeProject/
│
├── Final Home Proj.py      # Main party central! 🎉
├── README.md               # You are here! 📍
└── requirements.txt        # Shopping list for pip 🛒
```

### 🎛️ Key Parameters You Can Tweak

```python
# In Final Home Proj.py
class Ecoli:
    ATP = 200               # Starting energy
    glucose_sensitivity = 0.5    # How good they are at finding food
    base_growth_rate = 0.12     # How fast they multiply
```

Want to make changes? Here's what different values do:
- Higher ATP = Longer-lasting bacteria
- Higher sensitivity = Better at finding food
- Higher growth rate = Faster population growth

Remember: With great power comes great responsibility! Don't make your bacteria too powerful, or they might take over your CPU! 😅

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
