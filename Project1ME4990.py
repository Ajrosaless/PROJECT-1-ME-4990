"""
Conduction Studio - 1D Steady-State Heat Transfer GUI
For Heat Transfer Course - Chapter 3
Bergman, Lavine, Incropera, DeWitt - 6th Edition
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
from PIL import Image, ImageTk

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# Material Properties Dictionary (W/m·K)
MATERIALS = {
    "Copper": 385,
    "Aluminum": 205,
    "Steel (Carbon)": 54,
    "Stainless Steel": 15,
    "Concrete": 1.4,
    "Brick": 0.72,
    "Fiberglass Insulation": 0.04,
    "Foam Insulation": 0.03,
    "Glass": 1.4,
    "Wood (Oak)": 0.17,
    "Air Gap": 0.026,
}

# Thermal Contact Resistance Dictionary (m²·K/W × 10⁴)
CONTACT_RESISTANCE = {
    "Perfect Contact": 0,
    "Steel/Steel - Vacuum (100 kN/m²)": 6.0,
    "Steel/Steel - Vacuum (10000 kN/m²)": 0.7,
    "Aluminum/Aluminum - Vacuum (100 kN/m²)": 1.5,
    "Aluminum/Aluminum - Vacuum (10000 kN/m²)": 0.2,
    "Steel/Steel - Air": 2.75,
    "Aluminum/Aluminum - Air": 2.75,
    "Aluminum/Aluminum - Helium": 1.05,
    "Silicon Chip/Aluminum": 0.3,
    "Aluminum/Aluminum - Indium Foil": 0.07,
    "Stainless/Stainless - Dow Corning": 0.04,
}


class ConductionProblem:
    """Base class for all conduction problems"""

    def __init__(self):
        self.layers = []
        self.T_inner = 100
        self.T_outer = 20
        self.contact_resistances = []


class PlaneWall(ConductionProblem):
    """Composite plane wall with contact resistance"""

    def __init__(self):
        super().__init__()
        self.area = 1.0

    def add_layer(self, material, thickness):
        k = MATERIALS.get(material, 1.0)
        self.layers.append({'material': material, 'thickness': thickness, 'k': k})

    def add_contact_resistance(self, contact_type):
        R_tc = CONTACT_RESISTANCE.get(contact_type, 0) * 1e-4
        self.contact_resistances.append(R_tc)

    def calculate_heat_flux(self):
        R_total = sum(layer['thickness'] / (layer['k'] * self.area) for layer in self.layers)
        R_total += sum(R_tc / self.area for R_tc in self.contact_resistances)

        delta_T = self.T_inner - self.T_outer
        Q = delta_T / R_total
        q_flux = Q / self.area

        return q_flux, Q, R_total

    def get_temperature_profile(self):
        q_flux, Q, R_total = self.calculate_heat_flux()

        positions = [0]
        temperatures = [self.T_inner]
        current_x = 0
        current_T = self.T_inner

        for i, layer in enumerate(self.layers):
            delta_T_cond = Q * layer['thickness'] / (layer['k'] * self.area)
            current_x += layer['thickness']
            current_T -= delta_T_cond

            positions.append(current_x)
            temperatures.append(current_T)

            if i < len(self.contact_resistances):
                delta_T_contact = Q * self.contact_resistances[i] / self.area
                current_T -= delta_T_contact
                temperatures.append(current_T)
                positions.append(current_x)

        return np.array(positions), np.array(temperatures)


class Cylinder(ConductionProblem):
    """Composite cylindrical system"""

    def __init__(self):
        super().__init__()
        self.length = 1.0

    def add_layer(self, material, r_inner, r_outer):
        k = MATERIALS.get(material, 1.0)
        self.layers.append({'material': material, 'r_inner': r_inner, 'r_outer': r_outer, 'k': k})

    def add_contact_resistance(self, contact_type):
        R_tc = CONTACT_RESISTANCE.get(contact_type, 0) * 1e-4
        self.contact_resistances.append(R_tc)

    def calculate_heat_flux(self):
        R_total = 0

        for layer in self.layers:
            R_cond = np.log(layer['r_outer'] / layer['r_inner']) / (2 * np.pi * layer['k'] * self.length)
            R_total += R_cond

        for i, R_tc in enumerate(self.contact_resistances):
            if i < len(self.layers):
                r_interface = self.layers[i]['r_outer']
                A_interface = 2 * np.pi * r_interface * self.length
                R_total += R_tc / A_interface

        delta_T = self.T_inner - self.T_outer
        Q = delta_T / R_total

        return Q, R_total

    def get_temperature_profile(self):
        Q, R_total = self.calculate_heat_flux()

        if len(self.layers) == 0:
            return np.array([0]), np.array([self.T_inner])

        radii = []
        temperatures = []
        current_T = self.T_inner

        for i, layer in enumerate(self.layers):
            r1 = layer['r_inner']
            r2 = layer['r_outer']
            k = layer['k']

            r_layer = np.linspace(r1, r2, 20)

            for r in r_layer:
                T_r = current_T - (Q / (2 * np.pi * k * self.length)) * np.log(r / r1)
                radii.append(r)
                temperatures.append(T_r)

            current_T = temperatures[-1]
            if i < len(self.contact_resistances):
                A_interface = 2 * np.pi * r2 * self.length
                delta_T_contact = Q * self.contact_resistances[i] / A_interface
                current_T -= delta_T_contact

        return np.array(radii), np.array(temperatures)


class Sphere(ConductionProblem):
    """Sphere with internal heat generation"""

    def __init__(self):
        super().__init__()
        self.r_inner = 0.0
        self.r_outer = 0.1
        self.q_gen = 0
        self.h_outer = 0
        self.T_inf = 25
        self.material = "Copper"
        self.k = MATERIALS["Copper"]

    def set_geometry(self, r_inner, r_outer, material):
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.material = material
        self.k = MATERIALS.get(material, 1.0)

    def set_heat_generation(self, q_gen):
        self.q_gen = q_gen

    def set_convection(self, h_outer, T_inf):
        self.h_outer = h_outer
        self.T_inf = T_inf

    def calculate_heat_flux(self):
        if self.r_inner == 0:
            volume = (4 / 3) * np.pi * self.r_outer ** 3
        else:
            volume = (4 / 3) * np.pi * (self.r_outer ** 3 - self.r_inner ** 3)

        Q_gen = self.q_gen * volume
        A_outer = 4 * np.pi * self.r_outer ** 2
        q_flux = Q_gen / A_outer

        if self.h_outer > 0:
            T_surface = self.T_inf + q_flux / self.h_outer
        else:
            T_surface = self.T_outer

        if self.r_inner == 0:
            T_max = T_surface + (self.q_gen * self.r_outer ** 2) / (6 * self.k)
        else:
            T_max = T_surface + (self.q_gen / (6 * self.k)) * (self.r_outer ** 2 - self.r_inner ** 2)

        return q_flux, Q_gen, T_max, T_surface

    def get_temperature_profile(self):
        q_flux, Q_gen, T_max, T_surface = self.calculate_heat_flux()

        if self.r_inner == 0:
            r_points = np.linspace(0, self.r_outer, 100)
        else:
            r_points = np.linspace(self.r_inner, self.r_outer, 100)

        temperatures = []
        for r in r_points:
            T_r = T_surface + (self.q_gen / (6 * self.k)) * (self.r_outer ** 2 - r ** 2)
            temperatures.append(T_r)

        return r_points, np.array(temperatures)


class ConductionStudioApp(tk.Tk):
    """Main GUI application"""

    def __init__(self):
        super().__init__()
        self.title("Conduction Studio - 1D Steady-State Heat Transfer")
        self.geometry("1400x800")

        self.current_problem = None
        self.problem_type = tk.StringVar(value="Plane Wall")
        self.create_widgets()

    def create_widgets(self):
        main_frame = ttk.Frame(self)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.input_frame = ttk.LabelFrame(main_frame, text="Input Parameters", padding=10)
        self.input_frame.grid(row=0, column=0, sticky="nsew", padx=5)

        self.viz_frame = ttk.LabelFrame(main_frame, text="Visualization", padding=10)
        self.viz_frame.grid(row=0, column=1, sticky="nsew", padx=5)

        self.results_frame = ttk.LabelFrame(main_frame, text="Results", padding=10)
        self.results_frame.grid(row=0, column=2, sticky="nsew", padx=5)

        main_frame.columnconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=2)
        main_frame.columnconfigure(2, weight=1)
        main_frame.rowconfigure(0, weight=1)

        self.build_input_panel()
        self.build_viz_panel()
        self.build_results_panel()

    def build_input_panel(self):
        row = 0

        # --- Display header image using PIL (Pillow) ---
        try:
            # Load an image (replace 'banner.png' with your actual image file)
            img = Image.open("banner.png")
            img = img.resize((320, 120))  # Resize to fit the sidebar width
            self.tk_img = ImageTk.PhotoImage(img)  # Keep reference to avoid garbage collection

            # Create image label
            lbl_img = ttk.Label(self.input_frame, image=self.tk_img)
            lbl_img.grid(row=row, column=0, columnspan=3, pady=10)
            row += 1
        except Exception as e:
            print("Could not load banner image:", e)

        # --- Continue with the normal input layout ---
        ttk.Label(self.input_frame, text="Geometry Type:", font=('Arial', 10, 'bold')).grid(row=row, column=0,
                                                                                            sticky="w", pady=5)
        row += 1

        geometry_combo = ttk.Combobox(
            self.input_frame,
            textvariable=self.problem_type,
            values=["Plane Wall", "Cylinder", "Sphere"],
            state="readonly",
            width=25
        )
        geometry_combo.grid(row=row, column=0, columnspan=2, pady=5)
        geometry_combo.bind('<<ComboboxSelected>>', self.on_geometry_change)
        row += 1

        self.layers_label = ttk.Label(self.input_frame, text="Number of Layers:", font=('Arial', 9))
        self.layers_label.grid(row=row, column=0, sticky="w", pady=5)
        row += 1

        self.num_layers = tk.IntVar(value=2)
        self.layer_slider = ttk.Scale(
            self.input_frame,
            from_=1,
            to=5,
            orient=tk.HORIZONTAL,
            variable=self.num_layers,
            command=self.on_layer_change
        )
        self.layer_slider.grid(row=row, column=0, columnspan=2, sticky="ew", pady=5)

        self.layer_label = ttk.Label(self.input_frame, text="2 layers")
        self.layer_label.grid(row=row, column=2, pady=5)
        row += 1

        self.layer_inputs_frame = ttk.Frame(self.input_frame)
        self.layer_inputs_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1
        # Sphere frame
        self.sphere_frame = ttk.LabelFrame(self.input_frame, text="Sphere Parameters", padding=5)
        self.sphere_row = row
        row += 1

        ttk.Label(self.sphere_frame, text="Material:").grid(row=0, column=0, sticky="w", pady=2)
        self.sphere_material = tk.StringVar(value="Copper")
        ttk.Combobox(self.sphere_frame, textvariable=self.sphere_material,
                     values=list(MATERIALS.keys()), state="readonly", width=18).grid(row=0, column=1, padx=5, pady=2)

        ttk.Label(self.sphere_frame, text="Inner Radius (m):").grid(row=1, column=0, sticky="w", pady=2)
        self.sphere_r_inner = ttk.Entry(self.sphere_frame, width=12)
        self.sphere_r_inner.insert(0, "0.0")
        self.sphere_r_inner.grid(row=1, column=1, sticky="w", padx=5, pady=2)

        ttk.Label(self.sphere_frame, text="Outer Radius (m):").grid(row=2, column=0, sticky="w", pady=2)
        self.sphere_r_outer = ttk.Entry(self.sphere_frame, width=12)
        self.sphere_r_outer.insert(0, "0.05")
        self.sphere_r_outer.grid(row=2, column=1, sticky="w", padx=5, pady=2)

        ttk.Label(self.sphere_frame, text="Heat Gen q̇ (W/m³):").grid(row=3, column=0, sticky="w", pady=2)
        self.sphere_qgen = ttk.Entry(self.sphere_frame, width=12)
        self.sphere_qgen.insert(0, "1e6")
        self.sphere_qgen.grid(row=3, column=1, sticky="w", padx=5, pady=2)

        ttk.Label(self.sphere_frame, text="Convection h:").grid(row=4, column=0, sticky="w", pady=2)
        self.sphere_h = ttk.Entry(self.sphere_frame, width=12)
        self.sphere_h.insert(0, "100")
        self.sphere_h.grid(row=4, column=1, sticky="w", padx=5, pady=2)

        ttk.Label(self.sphere_frame, text="Ambient T∞ (°C):").grid(row=5, column=0, sticky="w", pady=2)
        self.sphere_Tinf = ttk.Entry(self.sphere_frame, width=12)
        self.sphere_Tinf.insert(0, "25")
        self.sphere_Tinf.grid(row=5, column=1, sticky="w", padx=5, pady=2)

        ttk.Separator(self.input_frame, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        self.temp_label = ttk.Label(self.input_frame, text="Temperatures:", font=('Arial', 10, 'bold'))
        self.temp_label.grid(row=row, column=0, sticky="w", pady=5)
        row += 1

        ttk.Label(self.input_frame, text="T₁ (inner/hot):").grid(row=row, column=0, sticky="w")
        self.T1_entry = ttk.Entry(self.input_frame, width=15)
        self.T1_entry.insert(0, "100")
        self.T1_entry.grid(row=row, column=1, sticky="w", padx=5)
        self.T1_row = row
        row += 1

        ttk.Label(self.input_frame, text="T₂ (outer/cold):").grid(row=row, column=0, sticky="w")
        self.T2_entry = ttk.Entry(self.input_frame, width=15)
        self.T2_entry.insert(0, "20")
        self.T2_entry.grid(row=row, column=1, sticky="w", padx=5)
        self.T2_row = row
        row += 1

        ttk.Separator(self.input_frame, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        button_frame = ttk.Frame(self.input_frame)
        button_frame.grid(row=row, column=0, columnspan=3, pady=10)

        ttk.Button(button_frame, text="Calculate", command=self.calculate).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Reset", command=self.reset).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Show Steps", command=self.show_solution_steps).pack(side=tk.LEFT, padx=5)

        self.layer_widgets = []
        self.update_layer_inputs()

    def build_viz_panel(self):
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.viz_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.ax.text(0.5, 0.5, 'Click Calculate', ha='center', va='center', fontsize=14, transform=self.ax.transAxes)
        self.ax.grid(True, alpha=0.3)
        self.canvas.draw()

    def build_results_panel(self):
        self.results_text = tk.Text(self.results_frame, width=35, height=30, wrap=tk.WORD)
        self.results_text.pack(fill=tk.BOTH, expand=True)
        self.results_text.insert('1.0', 'Results will appear here.')
        self.results_text.config(state=tk.DISABLED)

    def on_geometry_change(self, event=None):
        geom_type = self.problem_type.get()

        if geom_type == "Sphere":
            self.layers_label.grid_remove()
            self.layer_slider.grid_remove()
            self.layer_label.grid_remove()
            self.layer_inputs_frame.grid_remove()
            self.sphere_frame.grid(row=self.sphere_row, column=0, columnspan=3, sticky="ew", pady=10)
            self.temp_label.grid_remove()
            for w in self.input_frame.grid_slaves(row=self.T1_row): w.grid_remove()
            for w in self.input_frame.grid_slaves(row=self.T2_row): w.grid_remove()
        else:
            self.layers_label.grid()
            self.layer_slider.grid()
            self.layer_label.grid()
            self.layer_inputs_frame.grid()
            self.sphere_frame.grid_remove()
            self.temp_label.grid()
            for w in self.input_frame.grid_slaves(row=self.T1_row): w.grid()
            for w in self.input_frame.grid_slaves(row=self.T2_row): w.grid()

        self.update_layer_inputs()

    def on_layer_change(self, value):
        num = int(float(value))
        self.layer_label.config(text=f"{num} layers")
        self.update_layer_inputs()

    def update_layer_inputs(self):
        for widget in self.layer_inputs_frame.winfo_children():
            widget.destroy()

        self.layer_widgets = []
        geom_type = self.problem_type.get()

        if geom_type == "Sphere":
            return

        num_layers = self.num_layers.get()

        for i in range(num_layers):
            frame = ttk.LabelFrame(self.layer_inputs_frame, text=f"Layer {i + 1}", padding=5)
            frame.pack(fill=tk.X, pady=5)

            ttk.Label(frame, text="Material:").grid(row=0, column=0, sticky="w")
            material_var = tk.StringVar(value="Steel (Carbon)")
            ttk.Combobox(frame, textvariable=material_var, values=list(MATERIALS.keys()),
                         state="readonly", width=18).grid(row=0, column=1, padx=5)

            if geom_type == "Plane Wall":
                ttk.Label(frame, text="Thickness (m):").grid(row=1, column=0, sticky="w")
                thickness_entry = ttk.Entry(frame, width=10)
                thickness_entry.insert(0, "0.05")
                thickness_entry.grid(row=1, column=1, sticky="w", padx=5)
                self.layer_widgets.append({'material': material_var, 'thickness': thickness_entry})

            elif geom_type == "Cylinder":
                ttk.Label(frame, text="r_inner (m):").grid(row=1, column=0, sticky="w")
                r_inner_entry = ttk.Entry(frame, width=10)
                r_inner_entry.insert(0, str(0.02 + i * 0.01))
                r_inner_entry.grid(row=1, column=1, sticky="w", padx=5)

                ttk.Label(frame, text="r_outer (m):").grid(row=2, column=0, sticky="w")
                r_outer_entry = ttk.Entry(frame, width=10)
                r_outer_entry.insert(0, str(0.02 + (i + 1) * 0.01))
                r_outer_entry.grid(row=2, column=1, sticky="w", padx=5)
                self.layer_widgets.append(
                    {'material': material_var, 'r_inner': r_inner_entry, 'r_outer': r_outer_entry})

            if i < num_layers - 1:
                ttk.Label(frame, text="Contact R:").grid(row=3, column=0, sticky="w")
                contact_var = tk.StringVar(value="Perfect Contact")
                ttk.Combobox(frame, textvariable=contact_var, values=list(CONTACT_RESISTANCE.keys()),
                             state="readonly", width=18).grid(row=3, column=1, padx=5)
                self.layer_widgets[-1]['contact'] = contact_var

    def calculate(self):
        try:
            geom_type = self.problem_type.get()

            if geom_type == "Plane Wall":
                problem = PlaneWall()
                problem.T_inner = float(self.T1_entry.get())
                problem.T_outer = float(self.T2_entry.get())

                for wd in self.layer_widgets:
                    problem.add_layer(wd['material'].get(), float(wd['thickness'].get()))
                    if 'contact' in wd:
                        problem.add_contact_resistance(wd['contact'].get())

                q_flux, Q, R_total = problem.calculate_heat_flux()
                x_data, T_data = problem.get_temperature_profile()

                self.ax.clear()
                self.ax.plot(x_data, T_data, 'b-', linewidth=2)
                self.ax.scatter(x_data, T_data, c='red', s=30)
                self.ax.set_xlabel('Position x (m)')
                self.ax.set_ylabel('Temperature (°C)')
                self.ax.set_title('Plane Wall Temperature Distribution')
                self.ax.grid(True, alpha=0.3)
                self.canvas.draw()

                self.display_results_plane(q_flux, Q, R_total, problem)
                self.current_problem = problem

            elif geom_type == "Cylinder":
                problem = Cylinder()
                problem.T_inner = float(self.T1_entry.get())
                problem.T_outer = float(self.T2_entry.get())

                for wd in self.layer_widgets:
                    problem.add_layer(wd['material'].get(), float(wd['r_inner'].get()), float(wd['r_outer'].get()))
                    if 'contact' in wd:
                        problem.add_contact_resistance(wd['contact'].get())

                Q, R_total = problem.calculate_heat_flux()
                r_data, T_data = problem.get_temperature_profile()

                self.ax.clear()
                self.ax.plot(r_data, T_data, 'b-', linewidth=2)
                self.ax.set_xlabel('Radius r (m)')
                self.ax.set_ylabel('Temperature (°C)')
                self.ax.set_title('Cylinder Temperature Distribution')
                self.ax.grid(True, alpha=0.3)
                self.canvas.draw()

                self.display_results_cylinder(Q, R_total, problem)
                self.current_problem = problem

            elif geom_type == "Sphere":
                problem = Sphere()
                problem.set_geometry(float(self.sphere_r_inner.get()), float(self.sphere_r_outer.get()),
                                     self.sphere_material.get())
                problem.set_heat_generation(float(self.sphere_qgen.get()))
                problem.set_convection(float(self.sphere_h.get()), float(self.sphere_Tinf.get()))

                q_flux, Q_gen, T_max, T_surface = problem.calculate_heat_flux()
                r_data, T_data = problem.get_temperature_profile()

                self.ax.clear()
                self.ax.plot(r_data, T_data, 'b-', linewidth=2)
                self.ax.scatter([r_data[0], r_data[-1]], [T_data[0], T_data[-1]], c='red', s=50)
                self.ax.set_xlabel('Radius r (m)')
                self.ax.set_ylabel('Temperature (°C)')
                self.ax.set_title('Sphere with Heat Generation')
                self.ax.grid(True, alpha=0.3)
                self.canvas.draw()

                self.display_results_sphere(q_flux, Q_gen, T_max, T_surface, problem)
                self.current_problem = problem

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def display_results_plane(self, q_flux, Q, R_total, problem):
        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete('1.0', tk.END)
        self.results_text.insert('end', f"═══ PLANE WALL ═══\n\n")
        self.results_text.insert('end', f"Heat Flux: {q_flux:.2f} W/m²\n")
        self.results_text.insert('end', f"Heat Rate: {Q:.2f} W\n")
        self.results_text.insert('end', f"R_total: {R_total:.6f} K/W\n\n")
        for i, layer in enumerate(problem.layers):
            self.results_text.insert('end', f"Layer {i + 1}: {layer['material']}\n")
            self.results_text.insert('end', f"  k = {layer['k']:.2f} W/m·K\n")
            self.results_text.insert('end', f"  L = {layer['thickness']:.4f} m\n\n")
        self.results_text.config(state=tk.DISABLED)

    def display_results_cylinder(self, Q, R_total, problem):
        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete('1.0', tk.END)
        self.results_text.insert('end', f"═══ CYLINDER ═══\n\n")
        self.results_text.insert('end', f"Heat Rate: {Q:.2f} W\n")
        self.results_text.insert('end', f"R_total: {R_total:.6f} K/W\n\n")
        for i, layer in enumerate(problem.layers):
            self.results_text.insert('end', f"Layer {i + 1}: {layer['material']}\n")
            self.results_text.insert('end', f"  k = {layer['k']:.2f} W/m·K\n")
            self.results_text.insert('end', f"  r_i = {layer['r_inner']:.4f} m\n")
            self.results_text.insert('end', f"  r_o = {layer['r_outer']:.4f} m\n\n")
        self.results_text.config(state=tk.DISABLED)

    def display_results_sphere(self, q_flux, Q_gen, T_max, T_surface, problem):
        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete('1.0', tk.END)
        self.results_text.insert('end', f"═══ SPHERE ═══\n\n")
        self.results_text.insert('end', f"Material: {problem.material}\n")
        self.results_text.insert('end', f"k = {problem.k:.2f} W/m·K\n")
        self.results_text.insert('end', f"r_outer = {problem.r_outer:.4f} m\n")
        self.results_text.insert('end', f"q̇ = {problem.q_gen:.2e} W/m³\n\n")
        self.results_text.insert('end', f"Q_gen = {Q_gen:.2f} W\n")
        self.results_text.insert('end', f"q″_surface = {q_flux:.2f} W/m²\n")
        self.results_text.insert('end', f"T_max = {T_max:.2f} °C\n")
        self.results_text.insert('end', f"T_surface = {T_surface:.2f} °C\n")
        self.results_text.insert('end', f"ΔT = {T_max - T_surface:.2f} °C\n")
        self.results_text.config(state=tk.DISABLED)

    def reset(self):
        self.T1_entry.delete(0, tk.END)
        self.T1_entry.insert(0, "100")
        self.T2_entry.delete(0, tk.END)
        self.T2_entry.insert(0, "20")
        self.num_layers.set(2)
        self.update_layer_inputs()
        self.ax.clear()
        self.ax.text(0.5, 0.5, 'Click Calculate', ha='center', va='center', fontsize=14, transform=self.ax.transAxes)
        self.ax.grid(True, alpha=0.3)
        self.canvas.draw()
        self.current_problem = None

    def show_solution_steps(self):
        if self.current_problem is None:
            messagebox.showinfo("Info", "Calculate a problem first!")
            return

        window = tk.Toplevel(self)
        window.title("Step-by-Step Solution")
        window.geometry("700x600")

        frame = ttk.Frame(window)
        frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        scrollbar = ttk.Scrollbar(frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        text = tk.Text(frame, wrap=tk.WORD, yscrollcommand=scrollbar.set, font=('Courier', 9))
        text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=text.yview)

        geom = self.problem_type.get()

        if geom == "Plane Wall":
            self.generate_plane_wall_solution(text)
        elif geom == "Cylinder":
            self.generate_cylinder_solution(text)
        elif geom == "Sphere":
            self.generate_sphere_solution(text)

        text.config(state=tk.DISABLED)

    def generate_plane_wall_solution(self, txt):
        p = self.current_problem
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "STEP-BY-STEP SOLUTION: COMPOSITE PLANE WALL\n")
        txt.insert('end', "=" * 70 + "\n\n")

        txt.insert('end', "PROBLEM SETUP:\n")
        txt.insert('end', f"  Number of layers: {len(p.layers)}\n")
        txt.insert('end', f"  Inner temperature: T₁ = {p.T_inner}°C\n")
        txt.insert('end', f"  Outer temperature: T₂ = {p.T_outer}°C\n")
        txt.insert('end', f"  Area: A = {p.area} m²\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 1: Calculate Individual Thermal Resistances\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Conduction resistance formula for plane wall:\n")
        txt.insert('end', "  R_cond = L / (k × A)\n")
        txt.insert('end', "  where: L = thickness (m)\n")
        txt.insert('end', "         k = thermal conductivity (W/m·K)\n")
        txt.insert('end', "         A = cross-sectional area (m²)\n\n")

        R_tot = 0
        for i, layer in enumerate(p.layers):
            txt.insert('end', f"Layer {i + 1}: {layer['material']}\n")
            txt.insert('end', f"  Given: k = {layer['k']:.2f} W/m·K\n")
            txt.insert('end', f"         L = {layer['thickness']:.4f} m\n")
            txt.insert('end', f"         A = {p.area:.2f} m²\n")
            R = layer['thickness'] / (layer['k'] * p.area)
            txt.insert('end', f"  Calculation:\n")
            txt.insert('end', f"    R_cond,{i + 1} = L/(k×A)\n")
            txt.insert('end', f"    R_cond,{i + 1} = {layer['thickness']:.4f} / ({layer['k']:.2f} × {p.area:.2f})\n")
            txt.insert('end', f"    R_cond,{i + 1} = {R:.6f} K/W\n\n")
            R_tot += R

            if i < len(p.contact_resistances):
                Rc_area = p.contact_resistances[i]
                Rc = Rc_area / p.area
                txt.insert('end', f"Contact resistance between Layer {i + 1} and Layer {i + 2}:\n")
                txt.insert('end', f"  Given: R″_tc = {Rc_area:.6f} m²·K/W\n")
                txt.insert('end', f"  Calculation:\n")
                txt.insert('end', f"    R_tc = R″_tc / A\n")
                txt.insert('end', f"    R_tc = {Rc_area:.6f} / {p.area:.2f}\n")
                txt.insert('end', f"    R_tc = {Rc:.6f} K/W\n\n")
                R_tot += Rc

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 2: Calculate Total Thermal Resistance\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "For series resistances:\n")
        txt.insert('end', "  R_total = ΣR_cond + ΣR_contact\n")
        txt.insert('end', "  R_total = ")
        resistance_terms = []
        for i in range(len(p.layers)):
            resistance_terms.append(f"R_cond,{i + 1}")
            if i < len(p.contact_resistances):
                resistance_terms.append(f"R_tc,{i + 1}")
        txt.insert('end', " + ".join(resistance_terms) + "\n")
        txt.insert('end', f"  R_total = {R_tot:.6f} K/W\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 3: Calculate Heat Transfer Rate\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Using thermal circuit analogy (analogous to Ohm's Law):\n")
        txt.insert('end', "  Q = ΔT / R_total\n")
        txt.insert('end', "  where: ΔT = T₁ - T₂ (temperature difference)\n")
        txt.insert('end', "         R_total = total thermal resistance\n\n")

        delta_T = p.T_inner - p.T_outer
        Q = delta_T / R_tot
        txt.insert('end', f"  Given: T₁ = {p.T_inner}°C\n")
        txt.insert('end', f"         T₂ = {p.T_outer}°C\n")
        txt.insert('end', f"  Calculation:\n")
        txt.insert('end', f"    ΔT = T₁ - T₂ = {p.T_inner} - {p.T_outer} = {delta_T} K\n")
        txt.insert('end', f"    Q = ΔT / R_total\n")
        txt.insert('end', f"    Q = {delta_T} / {R_tot:.6f}\n")
        txt.insert('end', f"    Q = {Q:.2f} W\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 4: Calculate Heat Flux\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Heat flux formula:\n")
        txt.insert('end', "  q″ = Q / A\n")
        txt.insert('end', "  where: Q = heat transfer rate (W)\n")
        txt.insert('end', "         A = cross-sectional area (m²)\n\n")

        q_flux = Q / p.area
        txt.insert('end', f"  Calculation:\n")
        txt.insert('end', f"    q″ = Q / A\n")
        txt.insert('end', f"    q″ = {Q:.2f} / {p.area:.2f}\n")
        txt.insert('end', f"    q″ = {q_flux:.2f} W/m²\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 5: Calculate Interface Temperatures\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Temperature drop across each resistance:\n")
        txt.insert('end', "  ΔT = Q × R\n\n")

        current_T = p.T_inner
        txt.insert('end', f"Starting at left surface: T = {current_T:.2f}°C\n\n")

        for i, layer in enumerate(p.layers):
            R_layer = layer['thickness'] / (layer['k'] * p.area)
            delta_T_layer = Q * R_layer

            txt.insert('end', f"Through Layer {i + 1}:\n")
            txt.insert('end', f"  ΔT = Q × R_cond,{i + 1}\n")
            txt.insert('end', f"  ΔT = {Q:.2f} × {R_layer:.6f}\n")
            txt.insert('end', f"  ΔT = {delta_T_layer:.2f} K\n")
            current_T -= delta_T_layer
            txt.insert('end', f"  Temperature after Layer {i + 1}: T = {current_T:.2f}°C\n")

            if i < len(p.contact_resistances):
                R_tc = p.contact_resistances[i] / p.area
                delta_T_contact = Q * R_tc
                txt.insert('end', f"\n  Through contact resistance (temperature jump):\n")
                txt.insert('end', f"  ΔT_contact = Q × R_tc\n")
                txt.insert('end', f"  ΔT_contact = {Q:.2f} × {R_tc:.6f}\n")
                txt.insert('end', f"  ΔT_contact = {delta_T_contact:.2f} K\n")
                current_T -= delta_T_contact
                txt.insert('end', f"  Temperature after contact: T = {current_T:.2f}°C\n")
            txt.insert('end', "\n")

        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "FINAL ANSWERS:\n")
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', f"  Heat flux (q″)         = {q_flux:.2f} W/m²\n")
        txt.insert('end', f"  Heat transfer rate (Q) = {Q:.2f} W\n")
        txt.insert('end', f"  Total resistance (R)   = {R_tot:.6f} K/W\n")
        txt.insert('end', f"  Temperature drop (ΔT)  = {delta_T:.2f} K\n")
        txt.insert('end', "=" * 70 + "\n")

    def generate_cylinder_solution(self, txt):
        p = self.current_problem
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "STEP-BY-STEP SOLUTION: COMPOSITE CYLINDER\n")
        txt.insert('end', "=" * 70 + "\n\n")

        txt.insert('end', "PROBLEM SETUP:\n")
        txt.insert('end', f"  Number of layers: {len(p.layers)}\n")
        txt.insert('end', f"  Length: L = {p.length} m\n")
        txt.insert('end', f"  Inner temperature: T₁ = {p.T_inner}°C\n")
        txt.insert('end', f"  Outer temperature: T₂ = {p.T_outer}°C\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 1: Calculate Individual Thermal Resistances\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Conduction resistance formula for cylindrical geometry:\n")
        txt.insert('end', "  R_cond = ln(r₂/r₁) / (2πkL)\n")
        txt.insert('end', "  where: r₁ = inner radius (m)\n")
        txt.insert('end', "         r₂ = outer radius (m)\n")
        txt.insert('end', "         k = thermal conductivity (W/m·K)\n")
        txt.insert('end', "         L = cylinder length (m)\n\n")

        R_tot = 0.0
        R_layers_numeric = []  # keep numeric values for the summary line

        for i, layer in enumerate(p.layers):
            txt.insert('end', f"Layer {i + 1}: {layer['material']}\n")
            txt.insert('end', f"  Given: k = {layer['k']:.2f} W/m·K\n")
            txt.insert('end', f"         r₁ = {layer['r_inner']:.4f} m\n")
            txt.insert('end', f"         r₂ = {layer['r_outer']:.4f} m\n")
            txt.insert('end', f"         L = {p.length} m\n")

            ratio = layer['r_outer'] / layer['r_inner']
            ln_ratio = np.log(ratio)
            denom = 2 * np.pi * layer['k'] * p.length
            R = ln_ratio / denom

            txt.insert('end', "  Calculation:\n")
            txt.insert('end', f"    R_cond,{i + 1} = ln(r₂/r₁) / (2πkL)\n")
            txt.insert('end',
                       f"    R_cond,{i + 1} = ln({layer['r_outer']:.4f}/{layer['r_inner']:.4f}) / (2π × {layer['k']:.2f} × {p.length})\n")
            txt.insert('end', f"    R_cond,{i + 1} = ln({ratio:.4f}) / {denom:.4f}\n")
            txt.insert('end', f"    R_cond,{i + 1} = {ln_ratio:.6f} / {denom:.4f}\n")
            txt.insert('end', f"    R_cond,{i + 1} = {R:.6f} K/W\n\n")

            R_tot += R
            R_layers_numeric.append(R)

            # Contact resistance at outer radius of this layer (if provided)
            if i < len(p.contact_resistances):
                R_tc_area = p.contact_resistances[i]  # m²·K/W
                r_interface = layer['r_outer']  # m
                A_interface = 2 * np.pi * r_interface * p.length  # m²
                R_tc = R_tc_area / A_interface  # K/W

                txt.insert('end', f"Contact resistance between Layer {i + 1} and Layer {i + 2}:\n")
                txt.insert('end', f"  Given: R″_tc = {R_tc_area:.6f} m²·K/W\n")
                txt.insert('end', f"  Calculation:\n")
                txt.insert('end', f"    A_interface = 2πrL = 2π × {r_interface:.4f} × {p.length}\n")
                txt.insert('end', f"    A_interface = {A_interface:.4f} m²\n")
                txt.insert('end', f"    R_tc = R″_tc / A_interface = {R_tc_area:.6f} / {A_interface:.4f}\n")
                txt.insert('end', f"    R_tc = {R_tc:.6f} K/W\n\n")

                R_tot += R_tc
                R_layers_numeric.append(R_tc)

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 2: Calculate Total Thermal Resistance\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Total resistance formula for series resistances:\n")
        txt.insert('end', "  R_total = ΣR_cond + ΣR_contact\n\n")

        # Symbolic terms line (R_cond, R_tc …)
        terms = []
        for i in range(len(p.layers)):
            terms.append(f"R_cond,{i + 1}")
            if i < len(p.contact_resistances):
                terms.append(f"R_tc,{i + 1}")
        txt.insert('end', "  R_total = " + " + ".join(terms) + "\n")

        # Numeric values line (now using cylindrical formulas — no thickness/area!)
        txt.insert('end', "  R_total = " + " + ".join(f"{val:.6f}" for val in R_layers_numeric) + "\n")
        txt.insert('end', f"  R_total = {R_tot:.6f} K/W\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 3: Calculate Heat Transfer Rate\n")
        txt.insert('end', "-" * 70 + "\n\n")

        delta_T = p.T_inner - p.T_outer
        Q = delta_T / R_tot
        txt.insert('end', "Using thermal circuit analogy (Ohm’s law): Q = ΔT / R_total\n")
        txt.insert('end', f"  ΔT = {p.T_inner} - {p.T_outer} = {delta_T} K\n")
        txt.insert('end', f"  Q = {delta_T} / {R_tot:.6f} = {Q:.2f} W\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 4: Calculate Temperature Distribution\n")
        txt.insert('end', "-" * 70 + "\n\n")
        txt.insert('end', "Within each layer:  T(r) = T₁ - (Q/(2πkL)) ln(r/r₁)\n\n")

        current_T = p.T_inner
        for i, layer in enumerate(p.layers):
            R_layer = np.log(layer['r_outer'] / layer['r_inner']) / (2 * np.pi * layer['k'] * p.length)
            dT_layer = Q * R_layer
            T_inner_layer = current_T
            T_outer_layer = T_inner_layer - dT_layer

            txt.insert('end', f"Layer {i + 1} ({layer['r_inner']:.4f} m → {layer['r_outer']:.4f} m):\n")
            txt.insert('end', f"  ΔT_layer = Q × R_cond,{i + 1} = {Q:.2f} × {R_layer:.6f} = {dT_layer:.2f} K\n")
            txt.insert('end', f"  T(r₁) = {T_inner_layer:.2f}°C,  T(r₂) = {T_outer_layer:.2f}°C\n\n")

            current_T = T_outer_layer

            if i < len(p.contact_resistances):
                R_tc_area = p.contact_resistances[i]
                A_interface = 2 * np.pi * layer['r_outer'] * p.length
                R_tc = R_tc_area / A_interface
                dT_contact = Q * R_tc
                current_T -= dT_contact
                txt.insert('end',
                           f"Contact jump at r = {layer['r_outer']:.4f} m: ΔT = {dT_contact:.2f} K → T = {current_T:.2f}°C\n\n")

        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "FINAL ANSWERS:\n")
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', f"  Heat transfer rate (Q) = {Q:.2f} W\n")
        txt.insert('end', f"  Total resistance (R)   = {R_tot:.6f} K/W\n")
        txt.insert('end', f"  Temperature drop (ΔT)  = {delta_T:.2f} K\n")
        txt.insert('end', "=" * 70 + "\n")

    def generate_sphere_solution(self, txt):
        p = self.current_problem
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "STEP-BY-STEP SOLUTION: SPHERE WITH HEAT GENERATION\n")
        txt.insert('end', "=" * 70 + "\n\n")

        txt.insert('end', "PROBLEM SETUP:\n")
        if p.r_inner == 0:
            txt.insert('end', "  Type: Solid Sphere\n")
        else:
            txt.insert('end', "  Type: Hollow Sphere\n")
            txt.insert('end', f"  Inner radius: r_i = {p.r_inner:.4f} m\n")
        txt.insert('end', f"  Outer radius: r_o = {p.r_outer:.4f} m\n")
        txt.insert('end', f"  Material: {p.material}\n")
        txt.insert('end', f"  Thermal conductivity: k = {p.k:.2f} W/m·K\n")
        txt.insert('end', f"  Heat generation rate: q̇ = {p.q_gen:.2e} W/m³\n")
        txt.insert('end', f"  Convection coefficient: h = {p.h_outer:.2f} W/m²·K\n")
        txt.insert('end', f"  Ambient temperature: T∞ = {p.T_inf:.2f}°C\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 1: Calculate Total Heat Generation\n")
        txt.insert('end', "-" * 70 + "\n\n")

        if p.r_inner == 0:
            txt.insert('end', "Volume formula for solid sphere:\n")
            txt.insert('end', "  V = (4/3)πr³\n")
            txt.insert('end', "  where: r = sphere radius (m)\n\n")

            V = (4 / 3) * np.pi * p.r_outer ** 3
            txt.insert('end', f"  Given: r_o = {p.r_outer:.4f} m\n")
            txt.insert('end', f"  Calculation:\n")
            txt.insert('end', f"    V = (4/3)π({p.r_outer:.4f})³\n")
            txt.insert('end', f"    V = (4/3)π × {p.r_outer ** 3:.6f}\n")
            txt.insert('end', f"    V = {V:.6f} m³\n\n")
        else:
            txt.insert('end', "Volume formula for hollow sphere:\n")
            txt.insert('end', "  V = (4/3)π(r_o³ - r_i³)\n")
            txt.insert('end', "  where: r_o = outer radius (m)\n")
            txt.insert('end', "         r_i = inner radius (m)\n\n")

            V = (4 / 3) * np.pi * (p.r_outer ** 3 - p.r_inner ** 3)
            txt.insert('end', f"  Given: r_o = {p.r_outer:.4f} m\n")
            txt.insert('end', f"         r_i = {p.r_inner:.4f} m\n")
            txt.insert('end', f"  Calculation:\n")
            txt.insert('end', f"    V = (4/3)π({p.r_outer:.4f}³ - {p.r_inner:.4f}³)\n")
            txt.insert('end', f"    V = (4/3)π({p.r_outer ** 3:.6f} - {p.r_inner ** 3:.6f})\n")
            txt.insert('end', f"    V = (4/3)π × {p.r_outer ** 3 - p.r_inner ** 3:.6f}\n")
            txt.insert('end', f"    V = {V:.6f} m³\n\n")

        txt.insert('end', "Total heat generation:\n")
        txt.insert('end', "  Q_gen = q̇ × V\n")
        txt.insert('end', "  where: q̇ = volumetric heat generation rate (W/m³)\n")
        txt.insert('end', "         V = volume (m³)\n\n")

        Q_gen = p.q_gen * V
        txt.insert('end', f"  Given: q̇ = {p.q_gen:.2e} W/m³\n")
        txt.insert('end', f"  Calculation:\n")
        txt.insert('end', f"    Q_gen = {p.q_gen:.2e} × {V:.6f}\n")
        txt.insert('end', f"    Q_gen = {Q_gen:.2f} W\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 2: Calculate Heat Flux at Outer Surface\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "All generated heat must leave through outer surface.\n")
        txt.insert('end', "Surface area formula:\n")
        txt.insert('end', "  A_surface = 4πr_o²\n")
        txt.insert('end', "  where: r_o = outer radius (m)\n\n")

        A = 4 * np.pi * p.r_outer ** 2
        txt.insert('end', f"  Given: r_o = {p.r_outer:.4f} m\n")
        txt.insert('end', f"  Calculation:\n")
        txt.insert('end', f"    A_surface = 4π({p.r_outer:.4f})²\n")
        txt.insert('end', f"    A_surface = 4π × {p.r_outer ** 2:.6f}\n")
        txt.insert('end', f"    A_surface = {A:.4f} m²\n\n")

        txt.insert('end', "Heat flux formula:\n")
        txt.insert('end', "  q″_surface = Q_gen / A_surface\n")
        txt.insert('end', "  where: Q_gen = total heat generation (W)\n")
        txt.insert('end', "         A_surface = surface area (m²)\n\n")

        q_flux = Q_gen / A
        txt.insert('end', f"  Calculation:\n")
        txt.insert('end', f"    q″_surface = {Q_gen:.2f} / {A:.4f}\n")
        txt.insert('end', f"    q″_surface = {q_flux:.2f} W/m²\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 3: Calculate Surface Temperature\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Using convection boundary condition:\n")
        txt.insert('end', "  q″ = h(T_surface - T∞)\n")
        txt.insert('end', "  where: h = convection coefficient (W/m²·K)\n")
        txt.insert('end', "         T_surface = surface temperature (°C)\n")
        txt.insert('end', "         T∞ = ambient temperature (°C)\n\n")

        txt.insert('end', "Solving for T_surface:\n")
        txt.insert('end', "  T_surface = T∞ + q″/h\n\n")

        T_surf = p.T_inf + q_flux / p.h_outer
        txt.insert('end', f"  Given: T∞ = {p.T_inf:.2f}°C\n")
        txt.insert('end', f"         h = {p.h_outer:.2f} W/m²·K\n")
        txt.insert('end', f"         q″ = {q_flux:.2f} W/m²\n")
        txt.insert('end', f"  Calculation:\n")
        txt.insert('end', f"    T_surface = {p.T_inf:.2f} + {q_flux:.2f}/{p.h_outer:.2f}\n")
        txt.insert('end', f"    T_surface = {p.T_inf:.2f} + {q_flux / p.h_outer:.2f}\n")
        txt.insert('end', f"    T_surface = {T_surf:.2f}°C\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 4: Calculate Maximum Temperature\n")
        txt.insert('end', "-" * 70 + "\n\n")

        if p.r_inner == 0:
            txt.insert('end', "For solid sphere, T_max occurs at center (r = 0).\n")
            txt.insert('end', "Temperature distribution formula:\n")
            txt.insert('end', "  T(r) = T_surface + (q̇r_o²)/(6k) × [1 - (r/r_o)²]\n\n")
            txt.insert('end', "At r = 0:\n")
            txt.insert('end', "  T_max = T_surface + (q̇r_o²)/(6k)\n")
            txt.insert('end', "  where: q̇ = heat generation rate (W/m³)\n")
            txt.insert('end', "         r_o = outer radius (m)\n")
            txt.insert('end', "         k = thermal conductivity (W/m·K)\n\n")

            term = (p.q_gen * p.r_outer ** 2) / (6 * p.k)
            T_max = T_surf + term

            txt.insert('end', f"  Given: q̇ = {p.q_gen:.2e} W/m³\n")
            txt.insert('end', f"         r_o = {p.r_outer:.4f} m\n")
            txt.insert('end', f"         k = {p.k:.2f} W/m·K\n")
            txt.insert('end', f"  Calculation:\n")
            txt.insert('end', f"    T_max = T_surface + (q̇r_o²)/(6k)\n")
            txt.insert('end', f"    T_max = {T_surf:.2f} + ({p.q_gen:.2e} × {p.r_outer:.4f}²)/(6 × {p.k:.2f})\n")
            txt.insert('end', f"    T_max = {T_surf:.2f} + ({p.q_gen * p.r_outer ** 2:.4f})/{6 * p.k:.2f}\n")
            txt.insert('end', f"    T_max = {T_surf:.2f} + {term:.2f}\n")
            txt.insert('end', f"    T_max = {T_max:.2f}°C\n\n")
        else:
            txt.insert('end', "For hollow sphere, T_max occurs at inner surface.\n")
            txt.insert('end', "Temperature formula:\n")
            txt.insert('end', "  T_max = T_surface + (q̇/(6k)) × (r_o² - r_i²)\n")
            txt.insert('end', "  where: q̇ = heat generation rate (W/m³)\n")
            txt.insert('end', "         k = thermal conductivity (W/m·K)\n")
            txt.insert('end', "         r_o = outer radius (m)\n")
            txt.insert('end', "         r_i = inner radius (m)\n\n")

            term = (p.q_gen / (6 * p.k)) * (p.r_outer ** 2 - p.r_inner ** 2)
            T_max = T_surf + term

            txt.insert('end', f"  Given: q̇ = {p.q_gen:.2e} W/m³\n")
            txt.insert('end', f"         k = {p.k:.2f} W/m·K\n")
            txt.insert('end', f"         r_o = {p.r_outer:.4f} m\n")
            txt.insert('end', f"         r_i = {p.r_inner:.4f} m\n")
            txt.insert('end', f"  Calculation:\n")
            txt.insert('end',
                       f"    T_max = {T_surf:.2f} + ({p.q_gen:.2e}/(6 × {p.k:.2f})) × ({p.r_outer:.4f}² - {p.r_inner:.4f}²)\n")
            txt.insert('end',
                       f"    T_max = {T_surf:.2f} + ({p.q_gen / (6 * p.k):.4e}) × ({p.r_outer ** 2:.6f} - {p.r_inner ** 2:.6f})\n")
            txt.insert('end',
                       f"    T_max = {T_surf:.2f} + ({p.q_gen / (6 * p.k):.4e}) × {p.r_outer ** 2 - p.r_inner ** 2:.6f}\n")
            txt.insert('end', f"    T_max = {T_surf:.2f} + {term:.2f}\n")
            txt.insert('end', f"    T_max = {T_max:.2f}°C\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 5: Temperature Distribution\n")
        txt.insert('end', "-" * 70 + "\n\n")

        if p.r_inner == 0:
            txt.insert('end', "General temperature distribution for solid sphere:\n")
            txt.insert('end', "  T(r) = T_surface + (q̇/(6k)) × (r_o² - r²)\n\n")
            txt.insert('end', "This gives a parabolic temperature profile with\n")
            txt.insert('end', "maximum at the center (r = 0).\n\n")
        else:
            txt.insert('end', "General temperature distribution for hollow sphere:\n")
            txt.insert('end', "  T(r) = T_surface + (q̇/(6k)) × (r_o² - r²)\n")
            txt.insert('end', f"  Valid for {p.r_inner:.4f} ≤ r ≤ {p.r_outer:.4f} m\n\n")

        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "FINAL ANSWERS:\n")
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', f"  Total heat generation (Q_gen) = {Q_gen:.2f} W\n")
        txt.insert('end', f"  Surface heat flux (q″)        = {q_flux:.2f} W/m²\n")
        txt.insert('end', f"  Surface temperature (T_s)     = {T_surf:.2f}°C\n")
        txt.insert('end', f"  Maximum temperature (T_max)   = {T_max:.2f}°C\n")
        txt.insert('end', f"  Temperature difference (ΔT)   = {T_max - T_surf:.2f} K\n")
        txt.insert('end', "=" * 70 + "\n")
        p = self.current_problem
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "STEP-BY-STEP SOLUTION: COMPOSITE CYLINDER\n")
        txt.insert('end', "=" * 70 + "\n\n")

        txt.insert('end', "PROBLEM SETUP:\n")
        txt.insert('end', f"  Number of layers: {len(p.layers)}\n")
        txt.insert('end', f"  Length: L = {p.length} m\n")
        txt.insert('end', f"  Inner temperature: T₁ = {p.T_inner}°C\n")
        txt.insert('end', f"  Outer temperature: T₂ = {p.T_outer}°C\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 1: Calculate Individual Thermal Resistances\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Conduction resistance formula for cylindrical geometry:\n")
        txt.insert('end', "  R_cond = ln(r₂/r₁) / (2πkL)\n")
        txt.insert('end', "  where: r₁ = inner radius (m)\n")
        txt.insert('end', "         r₂ = outer radius (m)\n")
        txt.insert('end', "         k = thermal conductivity (W/m·K)\n")
        txt.insert('end', "         L = cylinder length (m)\n\n")

        R_tot = 0
        for i, layer in enumerate(p.layers):
            txt.insert('end', f"Layer {i + 1}: {layer['material']}\n")
            txt.insert('end', f"  Given: k = {layer['k']:.2f} W/m·K\n")
            txt.insert('end', f"         r₁ = {layer['r_inner']:.4f} m\n")
            txt.insert('end', f"         r₂ = {layer['r_outer']:.4f} m\n")
            txt.insert('end', f"         L = {p.length} m\n")

            R = np.log(layer['r_outer'] / layer['r_inner']) / (2 * np.pi * layer['k'] * p.length)

            txt.insert('end', f"  Calculation:\n")
            txt.insert('end', f"    R_cond,{i + 1} = ln(r₂/r₁) / (2πkL)\n")
            txt.insert('end',
                       f"    R_cond,{i + 1} = ln({layer['r_outer']:.4f}/{layer['r_inner']:.4f}) / (2π × {layer['k']:.2f} × {p.length})\n")
            txt.insert('end',
                       f"    R_cond,{i + 1} = ln({layer['r_outer'] / layer['r_inner']:.4f}) / {2 * np.pi * layer['k'] * p.length:.4f}\n")
            txt.insert('end',
                       f"    R_cond,{i + 1} = {np.log(layer['r_outer'] / layer['r_inner']):.6f} / {2 * np.pi * layer['k'] * p.length:.4f}\n")
            txt.insert('end', f"    R_cond,{i + 1} = {R:.6f} K/W\n\n")
            R_tot += R

            if i < len(p.contact_resistances):
                R_tc_area = p.contact_resistances[i]
                r_interface = layer['r_outer']
                A_interface = 2 * np.pi * r_interface * p.length
                R_tc = R_tc_area / A_interface

                txt.insert('end', f"Contact resistance at r = {r_interface:.4f} m:\n")
                txt.insert('end', f"  Given: R″_tc = {R_tc_area:.6f} m²·K/W\n")
                txt.insert('end', f"  Calculation:\n")
                txt.insert('end', f"    A_interface = 2πrL\n")
                txt.insert('end', f"    A_interface = 2π × {r_interface:.4f} × {p.length}\n")
                txt.insert('end', f"    A_interface = {A_interface:.4f} m²\n")
                txt.insert('end', f"    R_tc = R″_tc / A_interface\n")
                txt.insert('end', f"    R_tc = {R_tc_area:.6f} / {A_interface:.4f}\n")
                txt.insert('end', f"    R_tc = {R_tc:.6f} K/W\n\n")
                R_tot += R_tc

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 2: Calculate Total Thermal Resistance\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "For series resistances:\n")
        txt.insert('end', f"  R_total = ΣR_cond + ΣR_contact\n")
        txt.insert('end', f"  R_total = {R_tot:.6f} K/W\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 3: Calculate Heat Transfer Rate\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Using thermal circuit analogy:\n")
        txt.insert('end', "  Q = ΔT / R_total\n\n")

        delta_T = p.T_inner - p.T_outer
        Q = delta_T / R_tot

        txt.insert('end', f"  Given: T₁ = {p.T_inner}°C, T₂ = {p.T_outer}°C\n")
        txt.insert('end', f"  Calculation:\n")
        txt.insert('end', f"    ΔT = T₁ - T₂ = {p.T_inner} - {p.T_outer} = {delta_T} K\n")
        txt.insert('end', f"    Q = ΔT / R_total\n")
        txt.insert('end', f"    Q = {delta_T} / {R_tot:.6f}\n")
        txt.insert('end', f"    Q = {Q:.2f} W\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 4: Temperature Distribution\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Temperature varies logarithmically through each layer:\n")
        txt.insert('end', "  T(r) = T₁ - (Q/(2πkL)) × ln(r/r₁)\n\n")

        current_T = p.T_inner
        for i, layer in enumerate(p.layers):
            T_inner_layer = current_T
            coeff = Q / (2 * np.pi * layer['k'] * p.length)
            delta_T_layer = coeff * np.log(layer['r_outer'] / layer['r_inner'])
            T_outer_layer = T_inner_layer - delta_T_layer

            txt.insert('end', f"Layer {i + 1} ({layer['r_inner']:.4f} m < r < {layer['r_outer']:.4f} m):\n")
            txt.insert('end', f"  At r = {layer['r_inner']:.4f} m: T = {T_inner_layer:.2f}°C\n")
            txt.insert('end', f"  At r = {layer['r_outer']:.4f} m: T = {T_outer_layer:.2f}°C\n")
            txt.insert('end', f"  ΔT across layer = {delta_T_layer:.2f} K\n\n")

            current_T = T_outer_layer

            if i < len(p.contact_resistances):
                R_tc_area = p.contact_resistances[i]
                A_interface = 2 * np.pi * layer['r_outer'] * p.length
                delta_T_contact = Q * R_tc_area / A_interface
                current_T -= delta_T_contact
                txt.insert('end', f"  Temperature jump at interface: ΔT = {delta_T_contact:.2f} K\n")
                txt.insert('end', f"  Temperature after contact: T = {current_T:.2f}°C\n\n")

        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "FINAL ANSWERS:\n")
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', f"  Heat transfer rate (Q) = {Q:.2f} W\n")
        txt.insert('end', f"  Total resistance (R)   = {R_tot:.6f} K/W\n")
        txt.insert('end', f"  Temperature drop (ΔT)  = {delta_T:.2f} K\n")
        txt.insert('end', "=" * 70 + "\n")

    def gen_sphere_solution(self, txt):
        p = self.current_problem
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "STEP-BY-STEP SOLUTION: SPHERE WITH HEAT GENERATION\n")
        txt.insert('end', "=" * 70 + "\n\n")

        txt.insert('end', "PROBLEM SETUP:\n")
        if p.r_inner == 0:
            txt.insert('end', "  Type: Solid Sphere\n")
        else:
            txt.insert('end', "  Type: Hollow Sphere\n")
            txt.insert('end', f"  Inner radius: r_i = {p.r_inner:.4f} m\n")
        txt.insert('end', f"  Outer radius: r_o = {p.r_outer:.4f} m\n")
        txt.insert('end', f"  Material: {p.material}\n")
        txt.insert('end', f"  Thermal conductivity: k = {p.k:.2f} W/m·K\n")
        txt.insert('end', f"  Heat generation rate: q̇ = {p.q_gen:.2e} W/m³\n")
        txt.insert('end', f"  Convection coefficient: h = {p.h_outer:.2f} W/m²·K\n")
        txt.insert('end', f"  Ambient temperature: T∞ = {p.T_inf:.2f}°C\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 1: Calculate Total Heat Generation\n")
        txt.insert('end', "-" * 70 + "\n\n")

        if p.r_inner == 0:
            txt.insert('end', "Volume formula for solid sphere:\n")
            txt.insert('end', "  V = (4/3)πr³\n\n")
            V = (4 / 3) * np.pi * p.r_outer ** 3
            txt.insert('end', f"  Calculation:\n")
            txt.insert('end', f"    V = (4/3)π({p.r_outer:.4f})³\n")
            txt.insert('end', f"    V = {V:.6f} m³\n\n")
        else:
            txt.insert('end', "Volume formula for hollow sphere:\n")
            txt.insert('end', "  V = (4/3)π(r_o³ - r_i³)\n\n")
            V = (4 / 3) * np.pi * (p.r_outer ** 3 - p.r_inner ** 3)
            txt.insert('end', f"  Calculation:\n")
            txt.insert('end', f"    V = (4/3)π({p.r_outer:.4f}³ - {p.r_inner:.4f}³)\n")
            txt.insert('end', f"    V = (4/3)π({p.r_outer ** 3:.6f} - {p.r_inner ** 3:.6f})\n")
            txt.insert('end', f"    V = {V:.6f} m³\n\n")

        txt.insert('end', "Total heat generation:\n")
        txt.insert('end', "  Q_gen = q̇ × V\n")
        Q_gen = p.q_gen * V
        txt.insert('end', f"  Q_gen = {p.q_gen:.2e} × {V:.6f}\n")
        txt.insert('end', f"  Q_gen = {Q_gen:.2f} W\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 2: Calculate Heat Flux at Outer Surface\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "All generated heat must leave through outer surface.\n")
        txt.insert('end', "Surface area formula:\n")
        txt.insert('end', "  A_surface = 4πr_o²\n\n")

        A = 4 * np.pi * p.r_outer ** 2
        txt.insert('end', f"  Calculation:\n")
        txt.insert('end', f"    A_surface = 4π({p.r_outer:.4f})²\n")
        txt.insert('end', f"    A_surface = {A:.4f} m²\n\n")

        txt.insert('end', "Heat flux formula:\n")
        txt.insert('end', "  q″_surface = Q_gen / A_surface\n")
        q_flux = Q_gen / A
        txt.insert('end', f"  q″_surface = {Q_gen:.2f} / {A:.4f}\n")
        txt.insert('end', f"  q″_surface = {q_flux:.2f} W/m²\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 3: Calculate Surface Temperature\n")
        txt.insert('end', "-" * 70 + "\n\n")

        txt.insert('end', "Using convection boundary condition:\n")
        txt.insert('end', "  q″ = h(T_surface - T∞)\n")
        txt.insert('end', "Solving for T_surface:\n")
        txt.insert('end', "  T_surface = T∞ + q″/h\n\n")

        T_surf = p.T_inf + q_flux / p.h_outer
        txt.insert('end', f"  Calculation:\n")
        txt.insert('end', f"    T_surface = {p.T_inf:.2f} + {q_flux:.2f}/{p.h_outer:.2f}\n")
        txt.insert('end', f"    T_surface = {p.T_inf:.2f} + {q_flux / p.h_outer:.2f}\n")
        txt.insert('end', f"    T_surface = {T_surf:.2f}°C\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 4: Calculate Maximum Temperature\n")
        txt.insert('end', "-" * 70 + "\n\n")

        if p.r_inner == 0:
            txt.insert('end', "For solid sphere, T_max occurs at center (r = 0).\n")
            txt.insert('end', "Temperature distribution formula:\n")
            txt.insert('end', "  T(r) = T_surface + (q̇r_o²)/(6k) × [1 - (r/r_o)²]\n")
            txt.insert('end', "At r = 0:\n")
            txt.insert('end', "  T_max = T_surface + (q̇r_o²)/(6k)\n\n")

            T_max = T_surf + (p.q_gen * p.r_outer ** 2) / (6 * p.k)

            txt.insert('end', f"  Calculation:\n")
            txt.insert('end', f"    T_max = T_surface + (q̇r_o²)/(6k)\n")
            txt.insert('end', f"    T_max = {T_surf:.2f} + ({p.q_gen:.2e} × {p.r_outer:.4f}²)/(6 × {p.k:.2f})\n")
            txt.insert('end', f"    T_max = {T_surf:.2f} + ({p.q_gen * p.r_outer ** 2:.4f})/{6 * p.k:.2f}\n")
            txt.insert('end', f"    T_max = {T_surf:.2f} + {(p.q_gen * p.r_outer ** 2) / (6 * p.k):.2f}\n")
            txt.insert('end', f"    T_max = {T_max:.2f}°C\n\n")
        else:
            txt.insert('end', "For hollow sphere, T_max occurs at inner surface.\n")
            txt.insert('end', "Temperature formula:\n")
            txt.insert('end', "  T_max = T_surface + (q̇/(6k)) × (r_o² - r_i²)\n\n")

            T_max = T_surf + (p.q_gen / (6 * p.k)) * (p.r_outer ** 2 - p.r_inner ** 2)

            txt.insert('end', f"  Calculation:\n")
            txt.insert('end',
                       f"    T_max = {T_surf:.2f} + ({p.q_gen:.2e}/(6 × {p.k:.2f})) × ({p.r_outer:.4f}² - {p.r_inner:.4f}²)\n")
            txt.insert('end',
                       f"    T_max = {T_surf:.2f} + ({p.q_gen / (6 * p.k):.4e}) × ({p.r_outer ** 2 - p.r_inner ** 2:.6f})\n")
            txt.insert('end',
                       f"    T_max = {T_surf:.2f} + {(p.q_gen / (6 * p.k)) * (p.r_outer ** 2 - p.r_inner ** 2):.2f}\n")
            txt.insert('end', f"    T_max = {T_max:.2f}°C\n\n")

        txt.insert('end', "-" * 70 + "\n")
        txt.insert('end', "STEP 5: Temperature Distribution Equation\n")
        txt.insert('end', "-" * 70 + "\n\n")

        if p.r_inner == 0:
            txt.insert('end', "General temperature distribution for solid sphere:\n")
            txt.insert('end', "  T(r) = T_surface + (q̇/(6k)) × (r_o² - r²)\n\n")
            txt.insert('end', "This gives a parabolic temperature profile with\n")
            txt.insert('end', "maximum at the center.\n")
        else:
            txt.insert('end', "General temperature distribution for hollow sphere:\n")
            txt.insert('end', "  T(r) = T_surface + (q̇/(6k)) × (r_o² - r²)\n")
            txt.insert('end', f"  Valid for {p.r_inner:.4f} ≤ r ≤ {p.r_outer:.4f} m\n")

        txt.insert('end', "\n" + "=" * 70 + "\n")
        txt.insert('end', "FINAL ANSWERS:\n")
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', f"  Total heat generation (Q_gen) = {Q_gen:.2f} W\n")
        txt.insert('end', f"  Surface heat flux (q″)        = {q_flux:.2f} W/m²\n")
        txt.insert('end', f"  Surface temperature (T_s)     = {T_surf:.2f}°C\n")
        txt.insert('end', f"  Maximum temperature (T_max)   = {T_max:.2f}°C\n")
        txt.insert('end', f"  Temperature difference (ΔT)   = {T_max - T_surf:.2f} K\n")
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "=" * 70 + "\n\n")

        txt.insert('end', "PROBLEM SETUP:\n")
        txt.insert('end', f"  Layers: {len(p.layers)}\n")
        txt.insert('end', f"  T₁ = {p.T_inner}°C, T₂ = {p.T_outer}°C\n\n")

        txt.insert('end', "STEP 1: Calculate Thermal Resistances\n")
        txt.insert('end', "  R_cond = L/(k×A)\n\n")

        R_tot = 0
        for i, layer in enumerate(p.layers):
            txt.insert('end', f"Layer {i + 1} ({layer['material']}):\n")
            txt.insert('end', f"  k = {layer['k']:.2f} W/m·K\n")
            txt.insert('end', f"  L = {layer['thickness']:.4f} m\n")
            R = layer['thickness'] / (layer['k'] * p.area)
            txt.insert('end', f"  R = {R:.6f} K/W\n\n")
            R_tot += R

            if i < len(p.contact_resistances):
                Rc = p.contact_resistances[i] / p.area
                txt.insert('end', f"Contact R = {Rc:.6f} K/W\n\n")
                R_tot += Rc

        txt.insert('end', f"STEP 2: Total Resistance\n")
        txt.insert('end', f"  R_total = {R_tot:.6f} K/W\n\n")

        txt.insert('end', f"STEP 3: Heat Transfer Rate\n")
        txt.insert('end', f"  Q = ΔT/R = {p.T_inner - p.T_outer}/{R_tot:.6f}\n")
        Q = (p.T_inner - p.T_outer) / R_tot
        txt.insert('end', f"  Q = {Q:.2f} W\n\n")

        txt.insert('end', f"STEP 4: Heat Flux\n")
        txt.insert('end', f"  q″ = Q/A = {Q:.2f}/{p.area} = {Q / p.area:.2f} W/m²\n\n")

        txt.insert('end', "=" * 70 + "\n")

    def gen_cylinder_solution(self, txt):
        p = self.current_problem
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "STEP-BY-STEP SOLUTION: COMPOSITE CYLINDER\n")
        txt.insert('end', "=" * 70 + "\n\n")

        txt.insert('end', "PROBLEM SETUP:\n")
        txt.insert('end', f"  Layers: {len(p.layers)}\n")
        txt.insert('end', f"  T₁ = {p.T_inner}°C, T₂ = {p.T_outer}°C\n\n")

        txt.insert('end', "STEP 1: Cylindrical Resistances\n")
        txt.insert('end', "  R_cond = ln(r₂/r₁)/(2πkL)\n\n")

        R_tot = 0
        for i, layer in enumerate(p.layers):
            txt.insert('end', f"Layer {i + 1} ({layer['material']}):\n")
            txt.insert('end', f"  k = {layer['k']:.2f} W/m·K\n")
            txt.insert('end', f"  r₁ = {layer['r_inner']:.4f} m\n")
            txt.insert('end', f"  r₂ = {layer['r_outer']:.4f} m\n")
            R = np.log(layer['r_outer'] / layer['r_inner']) / (2 * np.pi * layer['k'] * p.length)
            txt.insert('end', f"  R = {R:.6f} K/W\n\n")
            R_tot += R

        txt.insert('end', f"STEP 2: Total Resistance\n")
        txt.insert('end', f"  R_total = {R_tot:.6f} K/W\n\n")

        txt.insert('end', f"STEP 3: Heat Transfer\n")
        Q = (p.T_inner - p.T_outer) / R_tot
        txt.insert('end', f"  Q = ΔT/R = {Q:.2f} W\n\n")

        txt.insert('end', "=" * 70 + "\n")

    def gen_sphere_solution(self, txt):
        p = self.current_problem
        txt.insert('end', "=" * 70 + "\n")
        txt.insert('end', "STEP-BY-STEP SOLUTION: SPHERE WITH HEAT GENERATION\n")
        txt.insert('end', "=" * 70 + "\n\n")

        txt.insert('end', "PROBLEM SETUP:\n")
        txt.insert('end', f"  Material: {p.material}\n")
        txt.insert('end', f"  k = {p.k:.2f} W/m·K\n")
        txt.insert('end', f"  r_o = {p.r_outer:.4f} m\n")
        txt.insert('end', f"  q̇ = {p.q_gen:.2e} W/m³\n")
        txt.insert('end', f"  h = {p.h_outer:.2f} W/m²·K\n")
        txt.insert('end', f"  T∞ = {p.T_inf:.2f}°C\n\n")

        txt.insert('end', "STEP 1: Total Heat Generation\n")
        V = (4 / 3) * np.pi * p.r_outer ** 3 if p.r_inner == 0 else (4 / 3) * np.pi * (p.r_outer ** 3 - p.r_inner ** 3)
        txt.insert('end', f"  V = {V:.6f} m³\n")
        Q_gen = p.q_gen * V
        txt.insert('end', f"  Q_gen = q̇×V = {Q_gen:.2f} W\n\n")

        txt.insert('end', "STEP 2: Surface Heat Flux\n")
        A = 4 * np.pi * p.r_outer ** 2
        q_flux = Q_gen / A
        txt.insert('end', f"  A = 4πr² = {A:.4f} m²\n")
        txt.insert('end', f"  q″ = Q/A = {q_flux:.2f} W/m²\n\n")

        txt.insert('end', "STEP 3: Surface Temperature\n")
        T_surf = p.T_inf + q_flux / p.h_outer
        txt.insert('end', f"  T_s = T∞ + q″/h\n")
        txt.insert('end', f"  T_s = {p.T_inf:.2f} + {q_flux:.2f}/{p.h_outer:.2f}\n")
        txt.insert('end', f"  T_s = {T_surf:.2f}°C\n\n")

        txt.insert('end', "STEP 4: Maximum Temperature\n")
        if p.r_inner == 0:
            txt.insert('end', "  T_max = T_s + (q̇r²)/(6k)\n")
            T_max = T_surf + (p.q_gen * p.r_outer ** 2) / (6 * p.k)
        else:
            txt.insert('end', "  T_max = T_s + (q̇/(6k))(r_o²-r_i²)\n")
            T_max = T_surf + (p.q_gen / (6 * p.k)) * (p.r_outer ** 2 - p.r_inner ** 2)
        txt.insert('end', f"  T_max = {T_max:.2f}°C\n\n")

        txt.insert('end', f"FINAL RESULTS:\n")
        txt.insert('end', f"  Q_gen = {Q_gen:.2f} W\n")
        txt.insert('end', f"  T_max = {T_max:.2f}°C\n")
        txt.insert('end', f"  T_surface = {T_surf:.2f}°C\n")
        txt.insert('end', f"  ΔT = {T_max - T_surf:.2f}°C\n")
        txt.insert('end', "=" * 70 + "\n\n")

        txt.insert('end', "Key Equations Summary:\n")
        txt.insert('end', "  • Volume:       V = (4/3)πr³  (solid) or (4/3)π(r_o³-r_i³)  (hollow)\n")
        txt.insert('end', "  • Heat gen:     Q = q̇ × V\n")
        txt.insert('end', "  • Surface area: A = 4πr_o²\n")
        txt.insert('end', "  • Heat flux:    q″ = Q/A\n")
        txt.insert('end', "  • Convection:   q″ = h(T_s - T∞)\n")
        txt.insert('end', "  • Max temp:     T_max = T_s + (q̇r_o²)/(6k)  (solid)\n")
        txt.insert('end', "=" * 70 + "\n")


if __name__ == "__main__":
    app = ConductionStudioApp()
    app.mainloop()

