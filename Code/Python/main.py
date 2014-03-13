from algos import *
from sorting import Sort as s
import numpy as np
import matplotlib.pyplot as plt
from runner import Runner
from grapher import TreeGrapher
from plotter import *

def tracefunc(frame, event, arg, indent=[0]):
  if event == "call":
    indent[0] += 2
    print "-" * indent[0] + "> call function", frame.f_code.co_name
  elif event == "return":
    print "<" + "-" * indent[0], "exit function", frame.f_code.co_name
    indent[0] -= 2
  return tracefunc

import sys
import warnings
warnings.filterwarnings('error')
#sys.settrace(tracefunc)

def f( x):
  if x == 0:
    return 0
  fx = np.sqrt(x)*np.log(x/1.02)
  return fx

def f(x):
  if x < 3.2:
    return x**2
  else:
    return x

#f = lambda x: x
boundary = [0, 10]

a = Binev2013( f, s, *boundary)
a.iterate( 18)
TreeGrapher( [a.tree_hp], 'lol.pdf', unique = False).graph()
plotfunc( f, *boundary)
plottree( a.tree_hp)
plt.show()

"""
sys.exit()
runner = Runner(n=100)
runner.add( Greedy          ( f, s, *args))
#runner.add( Binev2004First  ( f, s, 0, 1, 1, t = 0.01))
g = runner.add( Binev2004Second ( f, s, *args))
runner.add( Binev2007       ( f, s, *args))

runner.run( plot=True)
t = TreeGrapher( [g.tree], 'fux.pdf')
t.graph()
"""
