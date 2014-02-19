from algos import Binev2004First, Greedy
from sorting import Sort as s
import numpy as np

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
#inst = Greedy(f, s, 0, 1.0, 1)
inst = Binev2004First(f, s, 0.01, 0, 1, 1)
inst.iterate()
