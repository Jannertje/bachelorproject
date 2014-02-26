import matplotlib.pyplot as plt
from plotter import plotlists

class Runner( object):
  def __init__( self, n = False):
    self.insts = []
    self.n = n


  def add( self, algo):
    self.insts.append( algo)

  def run( self, plot = False):
    for inst in self.insts:
      inst.iterate( self.n)

    if plot:
      plotlists( [inst.error_list for inst in self.insts], 
                 [inst.__class__.__name__ for inst in self.insts])
      plt.show()
