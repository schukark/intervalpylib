from AreaCalculator import AreaCalculator
from SymbolicEquationSolver import SymbolicEquationSolver
import numpy as np
import intervalpy as ival

class Solver:
    def __init__(self, name):
        raise NotImplementedError

    def solve(self, system: SymbolicEquationSolver, boundary, verbose=False):
        raise NotImplementedError
    
class KrawczykSolver(Solver):
    def __init__(self, name):
        self._name = "Krawczyk"
    
    def __lambda_calcul(self, L, coef=1):
        mid_matrix = L
        n = L.shape[0]

        for i in range(n):
            for j in range(n):
                if not isinstance(mid_matrix[i, j], ival.Interval):
                    continue
                mid_matrix[i, j] = mid_matrix[i, j].mid()
        
        mid_matrix = mid_matrix.astype(dtype=np.float64)

        if np.linalg.det(mid_matrix) == 0:
            for i in range(n):
                mid_matrix[i, i] += 1
        
        Lambda = np.linalg.inv(mid_matrix)
        return Lambda
    
    def solve(self, system: SymbolicEquationSolver, grid, x_ini, verbose=False):
        inside_boxes = []
        border_boxes = []

        for box in grid:
            f = lambda x: (system.f_num_lam)(box, x)
            df = lambda x: (system.df_num_lam)(box, x)
            result = self.__iterative_method_root_exictence_default_ia(f, df, x_ini, self.__krawczyk, 
                                                                    self.__lambda_calcul, verbose=verbose)

            if result == 'I':
                inside_boxes.append(box)
            elif result == 'B':
                border_boxes.append(box)
        
        return inside_boxes, border_boxes

    def __krawczyk(self, f_num, df_num, x, c, Lambda, verbose=False):
        Kr = c - np.squeeze(np.matmul(Lambda, f_num(c))) + np.dot((np.identity(len(c)) - np.matmul(Lambda, df_num(x))), (x - c))

        for i in range(len(Kr)):
            if isinstance(Kr[i], ival.Interval):
                continue
            
            tmp = Kr[i]
            Kr[i] = ival.Interval([tmp, tmp])
        return Kr
    
    def __intersec(self, a, b):
        if a[1] < b[0] or b[1] < a[0]:
                return None
        else:
            return ival.Interval([max(a[0], b[0]), min(a[1], b[1])])
    
    def __iterative_method_root_exictence_default_ia(self, f, df, x_ini, method, lambda_calcul, 
                                                    eps = 10e-6, max_iters=100, verbose=True):
        n = len(x_ini)
        num_iters = 0
        queue = [x_ini]
        while len(queue)!=0:
            x = queue.pop(0)
            tol = max([x_i.width() for x_i in x])
            if np.all(ival.Interval([0, 0]).isIn(f[i](x[i])) for i in range(n)):
                num_iters+=1
                if verbose:
                    print("iter", num_iters)
                    print("\tqueue", queue)
                    print("\tx", x)
                c = np.array([xi.mid() for xi in x])
                Lambda = lambda_calcul(df(x))
                x_res = method(f, df, x, c, Lambda, verbose)
                if verbose:
                    print("\tMethod res", x_res)
                x_intersec = np.array([self.__intersec(x_res[i], x[i]) for i in range(n)])
                if verbose:
                    print("\tIntersection", x_intersec)
                if np.all([x_res[i].isIn(x[i]) for i in range(n)]):
                    return "I"
                elif np.any([type(x_i) == type(None) for x_i in x_intersec]) :
                    return "O"
                elif tol > eps and num_iters<max_iters:
                    queue.append(x_intersec)
            else:
                return "O"
        return "B"
    
if __name__ == "__main__":
    solver = KrawczykSolver("2-RPR")
    drawer = AreaCalculator("2-RPR", [3, 15, 8])
    grid = AreaCalculator.make_grid2d(left_bottom=(-20, -20), right_top=(20, 20), N=128)
    
    x_num = np.array([ival.Interval([3, 15])] * 2)
    inside, border = solver.solve(SymbolicEquationSolver("2-RPR", [8]), grid, x_num)

    drawer.uni_plotter(inside, border, size=2, ini_box=[[-20, 20]] * 2, title="Krawczyk")