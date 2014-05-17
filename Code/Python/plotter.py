import matplotlib.pyplot as plt
import numpy as np
from tree import Tree

def plotfunc( f, a, b = False):
  if isinstance(a, list):
    a, b = a
  X = np.linspace( a, b, 1024)
  plt.plot(X, [f(x) for x in X])

def plottree( tree):
  def islambda( v):
    return isinstance(v, type(lambda: None)) and v.__name__ == '<lambda>'
  for leaf in tree.leaves():
    if 'poly' in leaf.extra_info():
      p = leaf.extra_info_index('poly')
      if islambda( p):
        plotfunc( p, *leaf.boundary())
      else:
        plotfunc( lambda x: np.polyval( p[::-1], x), *leaf.boundary())

def plotlists( lists, legends, yscale = 'log'):
  for l in lists:
    plt.plot( l)
  plt.legend( legends)
  plt.yscale( yscale)
