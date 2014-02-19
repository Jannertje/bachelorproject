from abc import ABCMeta, abstractmethod
from tree import Tree_1D as t
from plotter import plottree, plotfunc
import matplotlib.pyplot as plt

class Algo(object):
  __metaclass__ = ABCMeta

  def __init__(self, f, s, a = 0, b = 1):
    #function to approximate
    self.f = f

    #sorting class
    self.s = s

    self.tree = t( -1., a, b)

  """
  Return the error of this node
  """
  @abstractmethod
  def error(self, node, d):
    pass

  @abstractmethod
  def needsSubdivide( self, leaf):
    pass

  def plot( self):
    plotfunc( self.f, *self.tree.boundary())
    plottree( self.tree)
    plt.show()

  def iteration( self):
    sorter = self.s()

    for leaf in self.tree.leaves():
      if self.needsSubdivide( leaf):
        error, info = self.error( leaf)
        leaf.value( error)
        leaf.extra_info( dict( leaf.extra_info().items() + info.items()))
        sorter.add( leaf)

    bests = sorter.find()
    return bests

  def iterate( self):
    l = []
    for i in range(10):
      bests = self.iteration()
      self.plot()
      l.append( self.tree.sum_of_leaves())
      if len(bests) == 0:
        break
      for leaf in bests:
        leaf.subdivide()

    plt.plot( l)
    plt.show()

class Binev2004First( Algo):
  def __init__( self, f, s, t, a = 0, b = 1, d = 1):
    self.d = d
    self.t = t
    super(self.__class__, self).__init__(f, s, a, b)

    from error_functionals import TwoNorm
    self.errorClass = TwoNorm(self.f)

  def __e( self, node, d):
    if node.value() > - 1.:
      return node.value()

    e, info = self.errorClass.error( node.boundary(), d)
    node.value( e)
    node.extra_info( dict( node.extra_info().items() + info.items()))
    return e

  def __lambda( self, node, d):
    if 'lambda' in node.extra_info():
      return node.extra_info()['lambda']

    #else: see (4.2) of Binev 2004
    sibling_errors = [self.__e( n, d) for n in node.siblings()]
    l = self.__e( node, d)/sum( sibling_errors)
    node.extra_info_index( 'lambda', l)
    return l

  def __alpha( self, node, d):
    if 'alpha' in node.extra_info():
      return node.extra_info()['alpha']

    #else: see (4.6) of Binev 2004
    if not node.getParent():
      return node.extra_info_index( 'alpha', 0.)

    l = self.__lambda( node, d)
    a = l * (self.__alpha( node.getParent(), d) \
           + self.__delta( node.getParent(), d))
    return node.extra_info_index( 'alpha', a)

  def __delta( self, node, d):
    #see (4.3) of Binev 2004
    def d():
      children_errors = [self.__e( n, d) for n in node.forest]
      return self.__e( node, d) - sum( children_errors)

    if 'delta' in node.extra_info():
      return node.extra_info()['delta']

    #else: see (4.4) of Binev 2004
    _d = d()
    if _d >= self.t:
      return node.extra_info_index( 'delta', 0.)
    else:
      return node.extra_info_index( 'delta', self.t - _d)

  def needsSubdivide( self, leaf):
    return self.error( leaf)[0] > self.t

  def error( self, node, d = False):
    if d == False:
      d = self.d
    
    if 'etilde' in node.extra_info():
      return node.extra_info()['etilde'], {}

    return node.extra_info_index( 
            'etilde', self.__e( node, d) - self.__alpha( node, d)
           ), {} #no extra info

class Greedy(Algo):
  def __init__(self, f, s, a = 0, b = 1, d = 1):
    self.d = d
    super(self.__class__, self).__init__(f, s, a, b)

    from error_functionals import TwoNorm
    self.errorClass = TwoNorm(self.f)

  def needsSubdivide( self, leaf):
    return True

  def error(self, node, d = False):
    if d == False:
      d = self.d
    v = node.value(), node.extra_info()
    if v[0] > -1.:
      return v
    return self.errorClass.error(node.boundary(), d)
