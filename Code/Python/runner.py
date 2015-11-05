import matplotlib.pyplot as plt
from plotter import plotlists

class Runner( object):
  def __init__( self, n = False):
    self.insts = []
    self.n = n


  def add( self, algo):
    self.insts.append( algo)
    return algo

  def inst_plot_string( self, inst):
    string = inst.__class__.__name__
    if hasattr( inst, 'd'):
      string += " d="
      string += str(inst.d)
    return string

  def run( self, plot = False):
    for inst in self.insts:
      inst.iterate( self.n)

    if plot:
      plotlists( [[inst.dof_list, inst.error_list] for inst in self.insts], 
                 [self.inst_plot_string( inst) for inst in self.insts])
      plt.show()
