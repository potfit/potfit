from geodesiclm import geodesiclm
import numpy as np


def func(x, A):
    """ Rosenbrock function """
    return np.array( [ (1 - x[0]), A*(x[1] - x[0]**2)] )


x0 = np.array([ -1.0, 1.0])
xf, info = geodesiclm(func, x0, args = (1000,), full_output=1, print_level = 5, iaccel = 1, maxiters = 10000, artol = -1.0, xtol = -1, ftol = -1, avmax = 2.0)

print info

print xf
