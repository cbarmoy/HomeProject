# 🦠 Bacterial Growth Simulation

A sophisticated simulation of bacterial growth with fluid dynamics, metabolite production, and interactive visualization.

## 🚀 Prerequisites

- Python 3.8 or higher
- Git (for cloning the repository)
- A modern web browser (for the dashboard)

## 🎮 Installation & Setup

1. **Clone the Repository**
   ```bash
   git clone <your-repo-url>
   cd HomeProject
   ```

2. **Create and Activate Virtual Environment**
   
   Windows:
   ```bash
   python -m venv .venv
   .venv\Scripts\activate
   ```
   
   Mac/Linux:
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   ```

3. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
   ```

## 🎯 Running the Simulation

1. **Start the Simulation**
   ```bash
   python src/main.py
   ```

2. **Access the Dashboard**
   - Open your web browser
   - Navigate to: `http://127.0.0.1:8050`
   - The dashboard will automatically update in real-time

## 📊 Features

### Real-time Visualization
- Interactive dashboard showing:
  - Bacterial positions and movements
  - Glucose distribution
  - Metabolite concentrations
  - Population statistics
  - ATP levels

### Simulation Components
- E. coli bacteria with ATP-based metabolism
- Glucose particles with physics-based movement
- Metabolite production and toxicity effects
- Fluid dynamics affecting particle behavior
- Drain system for particle collection and recycling

### Data Export
- Automatic CSV export of simulation data
- Time series of all key metrics
- Particle positions and states

## 🔧 Configuration

You can modify simulation parameters in `src/main.py`:

```python
scenarios = {
    'default': {
        'n_ecoli': 10,          # Number of initial bacteria
        'n_glucose': 300,       # Number of glucose particles
        'space_size': 20e-6,    # Size of simulation space
        'max_steps': 250,       # Simulation duration
        'glucose_feed_interval': 20,  # How often to add new glucose
        'glucose_feed_amount': 50     # How much glucose to add
    }
}
```

## 🧪 Advanced Usage

### Multiple Scenarios
The simulation supports running multiple scenarios:
```python
python src/main.py --scenario high_density
```

Available scenarios:
- `default`: Standard simulation
- `high_density`: More bacteria and glucose
- Custom scenarios can be added in `src/main.py`

### Debugging
For development and debugging:
```python
python src/main.py --debug
```

## 📁 Project Structure

```
src/
├── models/          # Particle classes (bacteria, glucose, metabolites)
├── physics/         # Fluid dynamics and spatial calculations
├── simulation/      # Core simulation logic
├── interface/       # Dashboard and visualization
├── systems/         # Drain system and recycling
└── utils/          # Helper functions and data export
```

## 🛠️ Troubleshooting

### Common Issues

1. **Import Errors**
   - Make sure you're running from the project root
   - Verify virtual environment is activated
   - Check all dependencies are installed

2. **Dashboard Not Loading**
   - Verify port 8050 is available
   - Check browser console for errors
   - Ensure all Dash dependencies are installed

3. **Performance Issues**
   - Reduce number of particles
   - Lower the simulation steps
   - Check system resource usage

## 📚 Dependencies

Key packages:
- numpy: Scientific computing
- dash: Interactive dashboard
- plotly: Data visualization
- pandas: Data handling
- scipy: Scientific calculations

See `requirements.txt` for complete list and versions.

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## 📝 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

- Built with Python and Dash
- Inspired by bacterial growth models
- Uses Navier-Stokes equations for fluid dynamics 