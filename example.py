from intervalpylib import KrawczykSolver, AreaCalculator, SymbolicEquationSolver
import intervalpy as ival
import numpy as np

solver = KrawczykSolver("2-RPR")
drawer = AreaCalculator("2-RPR", [3, 15, 8])
grid = AreaCalculator.make_grid2d(left_bottom=(-20, -20), right_top=(20, 20), N=128)

x_num = np.array([ival.Interval([3, 15])] * 2)
inside, border = solver.solve(SymbolicEquationSolver("2-RPR", [8]), grid, x_num)

drawer.uni_plotter(inside, border, size=2, ini_box=[[-20, 20]] * 2, title="Krawczyk")