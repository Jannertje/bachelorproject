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
  if x < 3: return np.exp( x)
  else: return 0.01

#f = lambda x: np.sin(x)
boundary = [0, 5]

runner = Runner( n = 5)
#a = Binev2004Second( f, s, *boundary, d=2)
a = Binev2013( f, s, *boundary)
runner.add( a)
#runner.add( Binev2007( f, s, *boundary, d=2))
#runner.add( Binev2007( f, s, *boundary, d=3))
runner.run( plot = True)
a.plot()

sys.exit()
"""
args = [0.5, 1, 3]
runner = Runner(n=5)
runner.add( Greedy          ( f, s, *args))
#runner.add( Binev2004First  ( f, s, 0, 1, 1, t = 0.01))
#runner.add( Binev2004Second ( f, s, *args))
#runner.add( Binev2007       ( f, s, *args))
#runner.add( Binev2013       ( f, s, *boundary))

runner.run( plot=True)
"""
