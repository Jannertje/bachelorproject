from error_functionals import TwoNormOrth as no
from error_functionals import TwoNorm as n
from plotter import plotfunc
import matplotlib.pyplot as plt
import numpy as np

f = lambda x: x**3
node = [0, 10]
e = n( f)
e2 = no( f)
pol, coeffs = e.bestpoly( node, d=17)
pol2, coeffs2 = e2.bestpoly( node, d=17)

plotfunc( lambda x: pol(x) - pol2(x), *node)
plt.show()
