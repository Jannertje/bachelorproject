import matplotlib.pyplot as plt
import numpy as np

def plotfunc( f, a, b = False):
  if isinstance(a, list):
    a, b = a
  X = np.linspace( a, b, 1024)
  plt.plot(X, [f(x) for x in X])
