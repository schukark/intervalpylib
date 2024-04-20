import numpy as np
import sympy


class SymbolicEquationSolver:
    def __init__(self, name, params):
        if name == "2-RPR":
            self._f, self._u, self._v = self.__symbolic_2rpr_func(*params)
        elif name == "3-RPR":
            self._f, self._u, self._v = self.__symbolic_3rpr_func(*params)
        elif name == "custom":
            assert len(params) >= 3
            self._f, self._u, self._v = params[0], params[1], params[2]
        
        self._f_lam = None
        self._df_lam = None
    
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

    def __derive_matrix(self, f, x):
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
    
    @property
    def f(self):
        return self._f
    
    @property
    def u(self):
        return self._u
    
    @property
    def v(self):
        return self._v

    @property
    def f_num_lam(self):
        if self._f_lam is None:
            self._f_lam = sympy.lambdify([self.u, self.v], self.f)
        return self._f_lam
    
    @property
    def df_num_lam(self):
        if self._df_lam is None:
            df = self.__derive_matrix(self.f, self.v)
            self._df_lam = sympy.lambdify([self.u, self.v], df)
        return self._df_lam
    
if __name__ == "__main__":
    two_rpr = SymbolicEquationSolver("2-RPR", [8])
    print(two_rpr.f, two_rpr.u, two_rpr.v)
