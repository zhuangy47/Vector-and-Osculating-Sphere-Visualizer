import numpy as np
from sympy import *
from sympy.physics.vector import *
import tkinter as tk
import tkinter.ttk as ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

class GraphInfo:
    def __init__(self, r, poi):
        self.r = r
        self.ori_poi = poi
        self.poi_with_offset = poi
        self.graph_artist = []
        self.quivers_artist = []
        self.plane_artist = []
        self.osc_sphere_artist = []
        self.graph_radius = 5
        self.draw_curvature = False
        self.draw_normal_plane = False
        self.draw_binormal_plane = False
        self.draw_tangent_plane = False
        self.draw_osculating_sphere = False
        self.r_gen = 0
        self.t_vec_gen = 0
        self.n_vec_gen = 0
        self.b_vec_gen = 0
        self.curvature_gen = 0
        self.r_poi_curr = []
        self.t_vec_curr = []
        self.n_vec_curr = []
        self.b_vec_curr = []
        self.curvature_curr = 0

def check_input(input: str):
    try:
        to_ret = sympify(input)
    except(Exception):
        print("Could not parse expression, try again.")
        raise Exception("Invalid Input")
    return to_ret

def set_axes_equal(ax: plt.Axes, x_points: list, y_points: list, z_points: list, origin: list):
    radius = round(0.5 * max(abs(max(x_points) - min(x_points)), abs(max(y_points) - min(y_points)), abs(max(z_points) - min(z_points))))
    ax.set_xlim3d(int(round(origin[0] - radius)), int(round(origin[0] + radius)))
    ax.set_ylim3d(int(round(origin[1] - radius)), int(round(origin[1] + radius)))
    ax.set_zlim3d(int(round(origin[2] - radius)), int(round(origin[2] + radius)))


def update_aux_slider(x, ax, info):
    # retrieve frequency
    f = float(x)
    info.poi_with_offset = info.ori_poi + f
    # update data
    graph_vectors(C, t, ax, info)
    draw_all_planes(ax, info)
    draw_osculating_sphere(ax, info)

    # required to update canvas and attached toolbar!

def full_redraw(ax, info: GraphInfo):
    # C = ReferenceFrame('C')
    # t = Symbol('t')
    vector_slider.set(0)
    validiate_input(info)
    graph_r_func(C, t, ax, info, info.graph_radius)
    graph_vectors(C, t, ax, info)
    draw_all_planes(ax, info)
    draw_osculating_sphere(ax, info)

def validiate_input(info: GraphInfo):
    components = []
    valid_r = []
    components.append(x_component.get())
    components.append(y_component.get())
    components.append(z_component.get())
    user_point = point_input.get()
    try:
        valid_poi = float(user_point)
    except:
        print("Please enter a float for your point")
        return
    for comp in components:
        try:
            sympified_input = check_input(comp)
            print(type(sympified_input))
        except(Exception):
            print(f"{comp} is not a valid expression")
            return
        valid_r.append(sympified_input)
    r_at_poi = [0, 0, 0]
    for dim, comp in enumerate(valid_r):
        if(isinstance(comp, int) or isinstance(comp, float)):
            r_at_poi[dim] = comp
        else:
            r_at_poi[dim] = comp.subs(t, valid_poi)
    # print(valid_r)
    # print(valid_poi)
    info.r = valid_r
    info.ori_poi = valid_poi
    info.poi_with_offset = valid_poi

def graph_r_func(C: ReferenceFrame, t: Symbol, ax, info: GraphInfo, radius):
    r_x = info.r[0]
    r_y = info.r[1]
    r_z = info.r[2]
    steps = np.linspace(info.ori_poi - radius, info.ori_poi + radius, int(resolution_slider.get()) * radius)
    steps = np.append(steps, info.ori_poi)
    steps = np.sort(steps)
    x_points = []
    y_points = []
    z_points = []
    r_at_poi = [0, 0, 0]
    for val in steps:
        calced_val = 0
        if(isinstance(r_x, int) or isinstance(r_x, float)):
            calced_val = r_x
        else:
            calced_val = r_x.subs(t, val)
        x_points.append(calced_val)
        if (val == info.ori_poi):
            r_at_poi[0] = calced_val

        if(isinstance(r_y, int) or isinstance(r_y, float)):
            calced_val = r_y
        else:
            calced_val = r_y.subs(t, val)
        y_points.append(calced_val)
        if (val == info.ori_poi):
            r_at_poi[1] = calced_val

        if(isinstance(r_z, int) or isinstance(r_z, float)):
            calced_val = r_z
        else:
            calced_val = r_z.subs(t, val)
        z_points.append(calced_val)
        if (val == info.ori_poi):
            r_at_poi[2] = calced_val

    for q in info.graph_artist:
        q[0].remove()
    info.graph_artist.clear()

    info.graph_artist.append(ax.plot3D(x_points, y_points, z_points, "blue"))
    ax.set_box_aspect([1,1,1])
    
    set_axes_equal(ax, x_points, y_points, z_points, r_at_poi)
    canvas.draw()


def graph_vectors(C: ReferenceFrame, t: Symbol, ax, info: GraphInfo):
    r_x = info.r[0]
    r_y = info.r[1]
    r_z = info.r[2]

    
    r_t = info.r[0] * C.x + info.r[1] * C.y + info.r[2] * C.z
    # print(r_t)
    info.r_gen = r_t

    # Calculate Tangent Vector
    t_t = r_t.diff(t, C)
    t_norm = t_t.normalize()
    info.t_vec_gen = t_norm

    # Calculate Normal Vector
    n_t = t_t.diff(t, C)
    n_norm = n_t.normalize()
    info.n_vec_gen = n_norm

    # Calcuate Binormal Vector
    b_t = n_t.cross(t_t)
    b_norm = b_t.normalize()
    info.b_vec_gen = b_norm

    ddt_t_norm = t_norm.diff(t, C)
    info.curvature_gen = sqrt(ddt_t_norm.dot(ddt_t_norm) / t_t.dot(t_t))
    

    # Get each vector's components and compute at poi
    t_norm_coeffs = t_norm.to_matrix(C)
    n_norm_coeffs = n_norm.to_matrix(C)
    b_norm_coeffs = b_norm.to_matrix(C)
    t_vec_at_poi = []
    n_vec_at_poi = []
    b_vec_at_poi = []
    curvature_at_poi = []
    for i in list(range(3)):
        t_vec_at_poi.append(N(t_norm_coeffs[i].subs(t, info.poi_with_offset)))
        n_vec_at_poi.append(N(n_norm_coeffs[i].subs(t, info.poi_with_offset)))
        b_vec_at_poi.append(-(N(b_norm_coeffs[i].subs(t, info.poi_with_offset))))
    curvature_at_poi.append(N(info.curvature_gen.subs(t, info.poi_with_offset)))
    
    # Get the value of the curve R at POI
    r_at_poi = [0, 0, 0]
    for dim, comp in enumerate(info.r):
        if(isinstance(comp, int) or isinstance(comp, float)):
            r_at_poi[dim] = comp
        else:
            r_at_poi[dim] = comp.subs(t, info.poi_with_offset)
    # print(f"tangent at t = {info.poi_with_offset} {t_vec_at_poi}")
    # print(f"normal at t = {info.poi_with_offset} {n_vec_at_poi}")
    # print(f"binormal at t = {info.poi_with_offset} {b_vec_at_poi}")
    # print(f"r at t = {info.poi_with_offset} {r_at_poi}")
    # print(f"curve at t = {info.poi_with_offset} {curvature_at_poi}")

    
    # Remove old vectors
    for q in info.quivers_artist:
        q.remove()
    info.quivers_artist.clear()
    
    # Graph current vectors
    info.quivers_artist.append(ax.quiver(*r_at_poi, *t_vec_at_poi, color = "green", arrow_length_ratio=0.1))
    info.quivers_artist.append(ax.quiver(*r_at_poi, *n_vec_at_poi, color = 'red', arrow_length_ratio=0.1))
    info.quivers_artist.append(ax.quiver(*r_at_poi, *b_vec_at_poi, color = 'orange', arrow_length_ratio=0.1))
    info.r_poi_curr = r_at_poi
    info.t_vec_curr = t_vec_at_poi
    info.n_vec_curr = n_vec_at_poi
    info.b_vec_curr = b_vec_at_poi
    info.curvature_curr = curvature_at_poi[0]

    canvas.draw()

def calc_x_plane(vector, y, z):
    return (((( float(vector[1]) / (0 - float(vector[0]))) * (y - float(info.r_poi_curr[1])) ) + (float(vector[2]) / (0 - float(vector[0]))) * (z - float(info.r_poi_curr[2]))) + float(info.r_poi_curr[0])) 

def calc_y_plane(vector, x, z):
    return ((( ( float(vector[0]) / (0 - float(vector[1]))) * (x - float(info.r_poi_curr[0])) ) + (float(vector[2]) / (0 - float(vector[1]))) * (z - float(info.r_poi_curr[2]))) + float(info.r_poi_curr[1]))


def calc_z_plane(vector, x, y):
    return ((( ( float(vector[0]) / (0 - float(vector[2]))) * (x - float(info.r_poi_curr[0])) ) + (float(vector[1]) / (0 - float(vector[2]))) * (y - float(info.r_poi_curr[1])))  + float(info.r_poi_curr[2]))

def calc_xy_plane(vector, x):
    return (( ( float(vector[0]) / (0 - float(vector[1]) )) * (x - float(info.r_poi_curr[0])) )  + float(info.r_poi_curr[1]))

def draw_all_planes(ax, info: GraphInfo):
    for q in info.plane_artist:
        q.remove()
    info.plane_artist.clear()
    if (info.draw_tangent_plane):
        draw_plane(ax, info.t_vec_curr, info, "green")
    if (info.draw_normal_plane):
        draw_plane(ax, info.n_vec_curr, info, "red")
    if (info.draw_binormal_plane):
        draw_plane(ax, info.b_vec_curr, info, "orange")
    canvas.draw()
    

def draw_plane(ax, vector: list[float], info: GraphInfo, input_color: str):
    

    if (vector[2] == 0):
        
        if (vector[1] == 0):
            y = np.linspace(float(info.r_poi_curr[1]) - 1, float(info.r_poi_curr[1]) + 1, 2)
            z = np.linspace(float(info.r_poi_curr[2]) - 1, float(info.r_poi_curr[2]) + 1, 2)
            Y, Z = np.meshgrid(y, z)
            X = Y*0 + float(info.r_poi_curr[0])
        elif(vector[0] == 0):
            x = np.linspace(float(info.r_poi_curr[0]) - 1, float(info.r_poi_curr[0]) + 1, 2)
            z = np.linspace(float(info.r_poi_curr[2]) - 1, float(info.r_poi_curr[2]) + 1, 2)
            X, Z = np.meshgrid(x, z)
            Y = X*0 + float(info.r_poi_curr[1])
        else:
            x = np.linspace(float(info.r_poi_curr[0]) - 1, float(info.r_poi_curr[0]) + 1, 2)
            z = np.linspace(float(info.r_poi_curr[2]) - 1, float(info.r_poi_curr[2]) + 1, 2)
            X, Z = np.meshgrid(x, z)
            Y = np.array(calc_xy_plane(vector, x))

    else:
        x = np.linspace(float(info.r_poi_curr[0]) - 1, float(info.r_poi_curr[0]) + 1, 2)
        y = np.linspace(float(info.r_poi_curr[1]) - 1, float(info.r_poi_curr[1]) + 1, 2)
        X, Y = np.meshgrid(x, y)

        Z = np.array(calc_z_plane(vector, X, Y))
    info.plane_artist.append(ax.plot_surface(X, Y, Z, color = input_color, alpha = 0.5))

def toggle_normal_plane(ax, info: GraphInfo):
    if normal_plane_toggle.config('relief')[-1] == 'sunken':
        normal_plane_toggle.config(relief="raised")
    else:
        normal_plane_toggle.config(relief="sunken")
    info.draw_normal_plane = not(info.draw_normal_plane)
    draw_all_planes(ax, info)

def toggle_binormal_plane(ax, info: GraphInfo):
    if binormal_plane_toggle.config('relief')[-1] == 'sunken':
        binormal_plane_toggle.config(relief="raised")
    else:
        binormal_plane_toggle.config(relief="sunken")
    info.draw_binormal_plane = not(info.draw_binormal_plane)
    draw_all_planes(ax, info)

def toggle_tangent_plane(ax, info: GraphInfo):
    if tangent_plane_toggle.config('relief')[-1] == 'sunken':
        tangent_plane_toggle.config(relief="raised")
    else:
        tangent_plane_toggle.config(relief="sunken")
    info.draw_tangent_plane = not(info.draw_tangent_plane)
    draw_all_planes(ax, info)

def set_graph_radius(val, ax, info):
    info.graph_radius = int(val)
    full_redraw(ax, info)

def toggle_osculating_sphere(ax, info: GraphInfo):
    if osculating_sphere_toggle.config('relief')[-1] == 'sunken':
        osculating_sphere_toggle.config(relief = 'raised')
    else:
        osculating_sphere_toggle.config(relief = 'sunken')
    info.draw_osculating_sphere = not(info.draw_osculating_sphere)
    
    draw_osculating_sphere(ax, info)
    canvas.draw()


def draw_osculating_sphere(ax, info: GraphInfo):
    
    for q in info.osc_sphere_artist:
        q.remove()
    info.osc_sphere_artist.clear()    
    if(info.draw_osculating_sphere):
        rho = float(1 / (info.curvature_curr))
        u, v = np.mgrid[0:2 * np.pi:15j, 0:np.pi:10j]
        x = rho * np.cos(u) * np.sin(v) + float((info.n_vec_curr[0]) * rho + info.r_poi_curr[0])
        y = rho * np.sin(u) * np.sin(v) + float((info.n_vec_curr[1]) * rho + info.r_poi_curr[1])
        z = rho * np.cos(v) + float(info.n_vec_curr[2] * rho  + info.r_poi_curr[2])
        info.osc_sphere_artist.append(ax.plot_surface(x, y, z, color="orange", alpha = 0.2, edgecolors='k', lw=0.2))


C = ReferenceFrame('C')
t = Symbol('t')
fig = Figure(figsize=(5, 4), dpi=100)
ax = fig.add_subplot(projection='3d')
ax.set_box_aspect([1,1,1])
set_axes_equal(ax, [0.5], [0.5], [0.5], [0.5, 0.5, 0.5])
ax.set_xlabel('$X$', fontsize=20)
ax.set_ylabel('$Y$', fontsize=20)

info = GraphInfo([0, 0, 0], 0)  # store r and poi

root = tk.Tk()
root.title("Visualizer")

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
component_insertion_frame = tk.Frame(root, height=100)
component_insertion_frame.grid_propagate(False)

graph_control = tk.Frame(root)

toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
toolbar.update()

button_quit = ttk.Button(master=root, text="Quit", command=root.destroy)
radius_slider = tk.Scale(graph_control, from_=1, to=20, orient=tk.HORIZONTAL, 
                            tickinterval=5, resolution=1, length = 100, label="Graph Radius", command = lambda val: set_graph_radius(val, ax, info))
radius_slider.set(5)
vector_slider = tk.Scale(graph_control, from_=-info.graph_radius, to=info.graph_radius, orient=tk.HORIZONTAL, 
                            tickinterval=1, resolution=0.1, length = 150, label="POI +-", command = lambda x: update_aux_slider(x, ax, info))
resolution_slider = tk.Scale(graph_control, from_=1, to=25, orient=tk.HORIZONTAL, 
                            tickinterval=5, resolution=1, length = 100, label="Graph Resolution", command = lambda x: graph_r_func(C, t, ax, info, info.graph_radius))
resolution_slider.set(10)



tangent_plane_toggle = tk.Button(graph_control, text="Toggle Tangent", command = lambda: toggle_tangent_plane(ax, info), relief = "raised")
normal_plane_toggle = tk.Button(graph_control, text="Toggle Normal", command = lambda: toggle_normal_plane(ax, info), relief = "raised")
binormal_plane_toggle = tk.Button(graph_control, text="Toggle Binormal", command = lambda: toggle_binormal_plane(ax, info), relief = "raised")
osculating_sphere_toggle = tk.Button(graph_control, text="Toggle Osc. Sphere", command = lambda: toggle_osculating_sphere(ax, info), relief = "raised")



x_component = ttk.Entry(component_insertion_frame, width=10)
x_component_label = ttk.Label(component_insertion_frame, text="x Component: ")
y_component = ttk.Entry(component_insertion_frame, width=10)
y_component_label = ttk.Label(component_insertion_frame, text="y Component: ")
z_component = ttk.Entry(component_insertion_frame, width=10)
z_component_label = ttk.Label(component_insertion_frame, text="z Component: ")

point_input = ttk.Entry(component_insertion_frame, width =10)
point_input_label = ttk.Label(component_insertion_frame, text = "POI: ")

submit_button = ttk.Button(component_insertion_frame, text = "Submit", command = lambda: full_redraw(ax, info))

button_quit.pack(side=tk.BOTTOM)


canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
x_component_label.pack(side=tk.LEFT, anchor='center')
x_component.pack(side=tk.LEFT, anchor='center')
y_component_label.pack(side=tk.LEFT, anchor='center')
y_component.pack(side=tk.LEFT, anchor='center')
z_component_label.pack(side=tk.LEFT, anchor='center')
z_component.pack(side=tk.LEFT, anchor='center')
point_input_label.pack(side=tk.LEFT, anchor='center')
point_input.pack(side=tk.LEFT, anchor='center')
submit_button.pack(side=tk.LEFT, anchor='center')

vector_slider.pack(side=tk.LEFT, anchor = 'center')
resolution_slider.pack(side=tk.LEFT, anchor = 'center')
radius_slider.pack(side=tk.LEFT, anchor = 'center')

tangent_plane_toggle.pack(side=tk.LEFT, anchor = 'center')
normal_plane_toggle.pack(side=tk.LEFT, anchor = 'center')
binormal_plane_toggle.pack(side=tk.LEFT, anchor = 'center')
osculating_sphere_toggle.pack(side=tk.LEFT, anchor = 'center')


toolbar.pack(side=tk.BOTTOM, fill=tk.X)

graph_control.pack(side = tk.BOTTOM)

component_insertion_frame.pack(side = tk.BOTTOM)
tk.mainloop()
