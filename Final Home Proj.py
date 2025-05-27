import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Ellipse, Circle
from matplotlib.widgets import Slider
import copy
import csv

# ====================
# --- Classes ---
# ====================

@dataclass
class Glucose:
    x: float
    y: float
    diameter: float = 0.0009e-6

    def diffuse(self, diffusion_strength, space_size):
        self.x += np.random.normal(0, diffusion_strength)
        self.y += np.random.normal(0, diffusion_strength)
        self.x, self.y = np.clip([self.x, self.y], 0, space_size)

    def __deepcopy__(self, memo):
        return Glucose(x=self.x, y=self.y, diameter=self.diameter)

@dataclass
class Ecoli:
    x: float
    y: float
    ATP: float = 100
    width: float = 2e-6
    height: float = 1e-6
    produced_metabolite: int = 0

    def move(self, glucose_list, step_size, noise, space_size, detection_radius):
        in_range = [g for g in glucose_list if np.hypot(self.x - g.x, self.y - g.y) <= detection_radius]
        if in_range:
            nearest = min(in_range, key=lambda g: (self.x - g.x) ** 2 + (self.y - g.y) ** 2)
            dx, dy = nearest.x - self.x, nearest.y - self.y
            norm = np.hypot(dx, dy)
            if norm > 0:
                dx /= norm
                dy /= norm
            ndx, ndy = np.random.uniform(-1, 1), np.random.uniform(-1, 1)
            self.x += step_size * ((1 - noise) * dx + noise * ndx)
            self.y += step_size * ((1 - noise) * dy + noise * ndy)
            self.x, self.y = np.clip([self.x, self.y], 0, space_size)

    def consume(self, glucose, atp_gain):
        threshold = (glucose.diameter / 2) + max(self.width, self.height) / 2
        if np.hypot(self.x - glucose.x, self.y - glucose.y) <= threshold:
            self.ATP += atp_gain
            return True
        return False

    def produce_metabolite(self, atp_cost):
        if self.ATP >= atp_cost:
            self.ATP -= atp_cost
            self.produced_metabolite += 1

    def can_divide(self, threshold):
        return self.ATP >= threshold

    def divide(self, space_size):
        self.ATP /= 2
        return Ecoli(
            x=np.clip(self.x + np.random.uniform(-1e-6, 1e-6), 0, space_size),
            y=np.clip(self.y + np.random.uniform(-1e-6, 1e-6), 0, space_size),
            ATP=self.ATP
        )

    def is_alive(self, threshold):
        return self.ATP >= threshold

    def __deepcopy__(self, memo):
        return Ecoli(
            x=self.x,
            y=self.y,
            ATP=self.ATP,
            width=self.width,
            height=self.height,
            produced_metabolite=self.produced_metabolite
        )

# ====================
# --- Gillespie algorithm ---
# ====================

def gillespie_step(ecoli, space_size):
    rates = {
        'divide': 0.1 if ecoli.can_divide(200) else 0,
        'die': 0.05 if not ecoli.is_alive(10) else 0,
        'idle': 1.0
    }
    total_rate = sum(rates.values())
    if total_rate == 0:
        return 'idle', None
    rand = np.random.uniform(0, total_rate)
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

def create_ecoli_and_glucose(n_ecoli, n_glucose, space_size):
    ecoli_list = [Ecoli(np.random.uniform(0, space_size), np.random.uniform(0, space_size))
                  for _ in range(n_ecoli)]
    glucose_list = [Glucose(np.random.uniform(0, space_size), np.random.uniform(0, space_size))
                    for _ in range(n_glucose)]
    return ecoli_list, glucose_list

def simulate(n_ecoli=10, n_glucose=300, space_size=20e-6, max_steps=250):
    ecoli_list, glucose_list = create_ecoli_and_glucose(n_ecoli, n_glucose, space_size)
    states = []

    for step in range(max_steps):
        for e in ecoli_list:
            e.move(glucose_list, step_size=1e-6, noise=0.5, space_size=space_size, detection_radius=1e-5)

        for g in glucose_list:
            g.diffuse(diffusion_strength=1e-7, space_size=space_size)

        new_glucose_list = []
        for g in glucose_list:
            if not any(e.consume(g, atp_gain=20) for e in ecoli_list):
                new_glucose_list.append(g)
        glucose_list = new_glucose_list

        for e in ecoli_list:
            e.produce_metabolite(atp_cost=10)

        new_bact = []
        survivors = []
        for e in ecoli_list:
            action, result = gillespie_step(e, space_size)
            if action == 'divide' and result:
                new_bact.append(result)
                survivors.append(e)
            elif action == 'idle':
                survivors.append(e)

        ecoli_list = survivors + new_bact
        states.append((copy.deepcopy(ecoli_list), copy.deepcopy(glucose_list)))

        if not ecoli_list:
            print(f"Arrêt étape {step + 1}: plus de bactéries.")
            break

    return states, space_size

# ====================
# --- Résumé graphique + CSV ---
# ====================

def plot_summary(states):
    steps = list(range(len(states)))
    tot_met = [sum(e.produced_metabolite for e in ecolis) for ecolis, _ in states]
    avg_atp = [np.mean([e.ATP for e in ecolis]) if ecolis else 0 for ecolis, _ in states]
    pop_size = [len(ecolis) for ecolis, _ in states]

    plt.figure(figsize=(10, 6))
    plt.plot(steps, tot_met, label="Métabolites produits", lw=2)
    plt.plot(steps, avg_atp, label="ATP moyen", lw=2)
    plt.plot(steps, pop_size, label="Population", lw=2)
    plt.xlabel("Étapes")
    plt.ylabel("Valeur")
    plt.title("Résumé de la simulation")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("recapitulatif_simulation.png")
    plt.close()

    with open("donnees_simulation.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Étape", "Métabolites produits", "ATP moyen", "Population"])
        for step, met, atp, pop in zip(steps, tot_met, avg_atp, pop_size):
            writer.writerow([step, met, atp, pop])

# ====================
# --- Animation ---
# ====================

def animate(states, space_size, interval=500):
    fig, (ax_sim, ax_ts) = plt.subplots(1, 2, figsize=(14, 6))
    plt.subplots_adjust(bottom=0.25)
    scale = 1e6
    space_um = space_size * scale

    xdata, y_met, y_atp, y_pop = [], [], [], []
    max_frames = len(states)

    # Courbes
    line_met, = ax_ts.plot([], [], lw=2, label='Métabolites')
    line_atp, = ax_ts.plot([], [], lw=2, label='ATP moyen')
    line_pop, = ax_ts.plot([], [], lw=2, label='Population')

    legend_handles = [
        Circle((0, 0), radius=0.1, edgecolor='blue', facecolor='cyan', alpha=0.7),
        Ellipse((0, 0), width=5, height=2, facecolor='red', edgecolor='black')
    ]

    def setup_axes():
        ax_sim.clear()
        ax_sim.set_xlim(0, space_um)
        ax_sim.set_ylim(0, space_um)
        ax_sim.set_xlabel('x')
        ax_sim.set_ylabel('y')
        ax_sim.set_title('Simulation spatiale (µm)')
        ax_sim.legend(legend_handles, ['Glucose', 'E. coli'], loc='upper right')

        ax_ts.set_xlim(0, max_frames)
        ax_ts.set_ylim(0, 1)  # temporaire
        ax_ts.set_ylabel('Valeur')
        ax_ts.set_xlabel('Étape')
        ax_ts.set_title('Évolution')

    setup_axes()

    # Slider interactif
    ax_slider = plt.axes([0.15, 0.1, 0.7, 0.03])
    slider = Slider(ax_slider, 'Progression', 0, max_frames - 1, valinit=0, valstep=1)

    def draw_frame(frame):
        ecolis, glucoses = states[frame]
        ax_sim.clear()
        setup_axes()

        for g in glucoses:
            ax_sim.add_patch(Circle((g.x * scale, g.y * scale), 0.1,
                                    edgecolor='blue', facecolor='cyan', alpha=0.7))
        for e in ecolis:
            c = plt.cm.Reds(min(e.ATP / 200, 1))
            ax_sim.add_patch(Ellipse((e.x * scale, e.y * scale), e.width * scale, e.height * scale,
                                     facecolor=c, edgecolor='black', lw=0.5))

        # Update graph data
        if frame >= len(xdata):
            xdata.append(frame)
            y_met.append(sum(e.produced_metabolite for e in ecolis))
            y_atp.append(np.mean([e.ATP for e in ecolis]) if ecolis else 0)
            y_pop.append(len(ecolis))

        line_met.set_data(xdata, y_met)
        line_atp.set_data(xdata, y_atp)
        line_pop.set_data(xdata, y_pop)

        max_y = max(max(y_met), max(y_atp), max(y_pop), 1)
        ax_ts.set_ylim(0, max_y * 1.1)
        ax_ts.legend()

        fig.canvas.draw_idle()

    def on_slider_change(val):
        draw_frame(int(slider.val))

    slider.on_changed(on_slider_change)

    # Animation automatique
    def update(frame):
        slider.set_val(frame)
        return []

    ani = FuncAnimation(fig, update, frames=max_frames, interval=interval, repeat=False)
    return ani

# ====================
# --- Main ---
# ====================

if __name__ == '__main__':
    states, L = simulate()
    ani = animate(states, L, interval=300)
    plot_summary(states)
    plt.show()
