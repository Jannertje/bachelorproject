from algos import *
from sorting import Sort as s
import numpy as np
import matplotlib.pyplot as plt
from runner import Runner

def tracefunc(frame, event, arg, indent=[0]):
  if event == "call":
    indent[0] += 2
    print "-" * indent[0] + "> call function", frame.f_code.co_name
  elif event == "return":
    print "<" + "-" * indent[0], "exit function", frame.f_code.co_name
    indent[0] -= 2
  return tracefunc

import sys
#sys.settrace(tracefunc)

f = lambda x: np.sqrt(x)*np.log(x/1.02)
args = [0.2, 5, 1]
runner = Runner(n=100)
runner.add( Greedy          ( f, s, *args))
#runner.add( Binev2004First  ( f, s, 0, 1, 1, t = 0.01))
runner.add( Binev2004Second ( f, s, *args))
runner.add( Binev2007       ( f, s, *args))

runner.run( plot=True)
