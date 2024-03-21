import numpy as np
import sympy


class SymbolicEquationSolver:
    def __init__(self, name, params):
        if name not in ["2-RPR", "3-RPR"]:
            raise NotImplementedError
        
        if name == "2-RPR":
            self._f, self._u, self._v = self.__symbolic_2rpr_func(*params)
        if name == "3-RPR":
            self._f, self._u, self._v = self.__symbolic_3rpr_func(*params)

    @property
    def f(self):
        return self._f
    
    @property
    def u(self):
        return self._u
    
    @property
    def v(self):
        return self._v
    
    def __symbolic_2rpr_func(self, d=8):
        """
        Creating symbol variables for 2-RPR system
        :return: symbolic eq. system,
                symbolic u (fixed boxes),
                symbolic v (checking boxes)
        """
        v = sympy.symbols("v1, v2")
        u = sympy.symbols("u1, u2")
        f = sympy.Matrix(
            [
                [v[0] ** 2 - (u[0] + 0.5 * d) ** 2 - u[1] ** 2],
                [v[1] ** 2 - (u[0] - 0.5 * d) ** 2 - u[1] ** 2],
            ]
        )
        return f, u, v

    def __symbolic_3rpr_func(self, x_c, y_c):
        """
        Creating symbol variables for 3-RPR system
        :return: symbolic eq. system,
                symbolic u (fixed boxes),
                symbolic v (checking boxes)
        """
        v = sympy.symbols("v1, v2, v3")
        u = sympy.symbols("u1, u2")
        f = sympy.Matrix(
            [
                [-v[0] ** 2 + (u[0] - x_c[0]) ** 2 + (u[1] - y_c[0]) ** 2],
                [-v[1] ** 2 + (u[0] - x_c[1]) ** 2 + (u[1] - y_c[1]) ** 2],
                [-v[2] ** 2 + (u[0] - x_c[2]) ** 2 + (u[1] - y_c[2]) ** 2],
            ]
        )
        return f, u, v

    def derive_matrix(f, x):
        """
        Function for calculating partial derivative of matrix g
        :param f : array to be derived
        :param x : variables for derivative
        :return gv: derived matrix
        """
        g_v_all = []
        for i in range(len(x)):
            g_v_all.append(sympy.diff(f, x[i]))  # Calculate derivative of G with respect to v
        gv = sympy.Matrix()
        for i in range(len(g_v_all)):
            gv = sympy.Matrix([gv, g_v_all[i]])
        gv = gv.reshape(f.shape[0], len(x)).T
        return gv
    
if __name__ == "__main__":
    two_rpr = SymbolicEquationSolver("2-RPR", [8])
    print(two_rpr.f, two_rpr.u, two_rpr.v)