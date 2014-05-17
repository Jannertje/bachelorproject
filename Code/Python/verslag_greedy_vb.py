from algos import *
from sorting import Sort as s
import numpy as np
import matplotlib.pyplot as plt
from runner import Runner
from grapher import TreeGrapher
from optimal import OptimalTree
from plotter import *

def f( x):
  if x < 1:
    return np.sin(16 * np.pi * x)
  elif x < 1.5:
    return 0
  else:
    return 0.5

boundary = [0,2]

plotfunc(f, *boundary)
plt.show()

"""
a = Greedy( f, s, *boundary, d=1)
a.iterate(5)
TreeGrapher( [a.tree], 'verslag_greedy_vb.pdf', unique = False).graph()
"""

t = OptimalTree( f, *boundary, d = 1, graph = 1)
trees = t.iterate( 5)
TreeGrapher( [trees[-1][0]], 'bluh.pdf', unique = False).graph()
plottree( trees[-1][0])
plotfunc( f, *boundary)
plt.show()
