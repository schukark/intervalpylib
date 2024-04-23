from intervalpylib import AreaCalculator, SymbolicEquationSolver, KrawczykSolver
import intervalpy as ival
import numpy as np

config = "3-RPR"

if config == "2-RPR":
    # Choose solver
    solver = KrawczykSolver("2-RPR")
    # Choose drawer and its parameters (depend on the configuration)
    drawer = AreaCalculator("2-RPR", [3, 15, 8])
    # Choose grid
    grid = AreaCalculator.make_grid2d(left_bottom=(-20, -20), right_top=(20, 20), N=128)

    # Choose lengths of rods in the 2-RPR configuration
    x_num = np.array([ival.Interval([3, 15])] * 2)
    # Solve the area
    inside, border = solver.solve(SymbolicEquationSolver("2-RPR", [8]), grid, x_num)

    # Plot the computed area
    drawer.uni_plotter(inside, border, size=2, ini_box=[[-20, 20]] * 2, title="Krawczyk")
elif config == "3-RPR":
    # Choose solver
    solver = KrawczykSolver("3-RPR")
    
    # Choose parameters (this code just rotates the circles around the point 150 degrees)
    import math
    phi = math.radians(150)
    x_a = [-15, 15, 0]
    y_a = [-5 * np.sqrt(3), -5 * np.sqrt(3), 10 * np.sqrt(3)]
    x_b = [-5, 5, 0]
    y_b = [-5 * np.sqrt(3) / 3, -5 * np.sqrt(3) / 3, 10 * np.sqrt(3) / 3]
    x_c = np.zeros(3)
    y_c = np.zeros(3)

    for i in range(3):
        x_c[i] = x_a[i] - x_b[i] * np.cos(phi) + y_b[i] * np.sin(phi)
        y_c[i] = y_a[i] - x_b[i] * np.sin(phi) - y_b[i] * np.cos(phi)
    
    # Choose drawer
    drawer = AreaCalculator("3-RPR", [x_c, y_c, 12, 27])
    # Choose grid
    grid = AreaCalculator.make_grid2d(left_bottom=(-20, -20), right_top=(20, 20), N=128)

    # Choose lengths of rods in the 2-RPR configuration
    x_num = np.array([ival.Interval([12, 27])] * 3)
    # Solve the area
    inside, border = solver.solve(SymbolicEquationSolver("3-RPR", [x_c, y_c]), grid, x_num)

    # Plot the computed area
    drawer.uni_plotter(inside, border, size=2, ini_box=[[-20, 20]] * 2, title="Krawczyk")