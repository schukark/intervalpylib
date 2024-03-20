import intervalpy as ival
from sympy import lambdify
from sympy.abc import x
from typing import Tuple, List

def f():
    return x ** 2 - 2


def bisect(interval: ival.Interval) -> Tuple[ival.Interval, ival.Interval]:
    mid_point = interval.mid()
    return (ival.Interval([interval[0], mid_point]), ival.Interval([mid_point, interval[1]]))

def split_into_2(interval1: ival.Interval, interval2: ival.Interval) -> List[ival.Interval]:
        result1 = None
        if interval1 is not None:
            result1: List[ival.Interval] = solve(interval1)
        
        result2 = None
        if interval2 is not None:
            result2: List[ival.Interval] = solve(interval2)

        if result1 is None:
            return result2
        if result2 is None:
            return result1
        
        return result1 + result2

def solve(search_space: ival.Interval) -> List[ival.Interval]:
    current_approx: ival.Interval = search_space
    eps: float = 10 ** (-7)

    func_lam = lambdify(x, f())
    func_lam_prime = lambdify(x, f().diff(x))

    while (func_lam(current_approx).width() > eps):
        new_approx = current_approx.mid() - func_lam(current_approx.mid()) / func_lam_prime(current_approx)
        not_intersecting = new_approx.isNoIntersec(current_approx)
        #print(current_approx, new_approx)

        if isinstance(new_approx, ival.ExtendedInterval):
            result1 = None
            if not not_intersecting[0]:
                result1 = new_approx[0].intersec(current_approx)
            
            result2 = None
            if not not_intersecting[1]:
                result2 = new_approx[1].intersec(current_approx)
            
            return split_into_2(result1, result2)
        elif isinstance(new_approx, ival.Interval):
            if not_intersecting:
                return None
            if new_approx.intersec(current_approx) == current_approx:
                intervals: Tuple[ival.Interval] = bisect(current_approx)
                return split_into_2(intervals[0], intervals[1])
            current_approx = new_approx.intersec(current_approx)
            continue

    return [current_approx]
    
if __name__ == "__main__":
    result = solve(ival.Interval([-2.5, 2.5]))
    print(result)