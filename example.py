from intervalpylib import AreaCalculator, SymbolicEquationSolver, KrawczykSolver
import intervalpy as ival
import numpy as np

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