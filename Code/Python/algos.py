from abc import ABCMeta, abstractmethod
from tree import Tree_1D as t
from tree import Tree_1D_hp
from plotter import plottree, plotfunc
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
from grapher import TreeGrapher

class Algo(object):
  __metaclass__ = ABCMeta

  def __init__(self, f, s, a = 0, b = 1):
    #function to approximate
    self.f = f

    #sorting class
    self.s = s

    self.tree = t( -1., a, b)

    self.error_list = []

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

  def iterate( self, n = False):
    def run():
      bests = self.iteration()
      #self.plot()
      current_error = self.tree.sum_of_leaves()
      self.error_list.append( current_error)
      if len(bests) == 0:
        return True
      for leaf in bests:
        leaf.subdivide()

      return False

    if n == False:
      while True:
        if run():
          break
    else:
      for i in range(n):
        if run():
          break

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

class Binev( Algo):
  def needsSubdivide( self, leaf):
    return True

  def _e( self, node, d):
    if node.value() > - 1.:
      return node.value()

    e, info = self.errorClass.error( node.boundary(), d)
    node.value( e)
    node.extra_info( dict( node.extra_info().items() + info.items()))
    return e

  @abstractmethod
  def error( self, node, d = False):
    pass

class Binev2004First( Binev):
  def __init__( self, f, s, a = 0, b = 1, d = 1, t = 0.01):
    self.d = d
    self.t = t
    super(self.__class__, self).__init__(f, s, a, b)

    from error_functionals import TwoNorm
    self.errorClass = TwoNorm(self.f)

  def __lambda( self, node, d):
    if 'lambda' in node.extra_info():
      return node.extra_info()['lambda']

    #else: see (4.2) of Binev 2004
    sibling_errors = [self._e( n, d) for n in node.siblings()]
    l = self._e( node, d)/sum( sibling_errors)
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
      children_errors = [self._e( n, d) for n in node.forest]
      return self._e( node, d) - sum( children_errors)

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
            'etilde', self._e( node, d) - self.__alpha( node, d)
           ), {} #no extra info

class Binev2004Second( Binev):
  def __init__( self, f, s, a=0, b=1, d=1):
    self.d = d
    super( self.__class__, self).__init__( f, s, a, b)

    from error_functionals import TwoNorm
    self.errorClass = TwoNorm( self.f)

  def __q( self, node, d):
    if 'q' in node.extra_info():
      return node.extra_info()['q']

    #else: see (5.3) of Binev 2004
    etilde, info = self.error( node, d)
    children_errors = [self._e( n, d) for n in node.forest]
    q = sum( children_errors)/( self._e( node, d) + etilde) * etilde
    node.extra_info_index( 'q', q)
    return q

  def error( self, node, d = False):
    if d == False:
      d = self.d
    
    if 'etilde' in node.extra_info():
      return node.extra_info()['etilde'], {}

    #see line above (5.2) of Binev 2004
    if not node.getParent():
      return node.extra_info_index( 'etilde', self._e( node, d)), {}

    #see (5.2) of Binev 2004
    return node.extra_info_index( 
            'etilde', self.__q( node.getParent(), d)
           ), {} #no extra info

class Binev2007( Binev):
  def __init__( self, f, s, a=0, b=1, d=1):
    self.d = d
    super( self.__class__, self).__init__( f, s, a, b)

    from error_functionals import TwoNorm
    self.errorClass = TwoNorm( self.f)

  def error( self, node, d = False):
    if d == False:
      d = self.d
    
    if 'etilde' in node.extra_info():
      return node.extra_info()['etilde'], {}

    #see (4) of Binev 2007
    if not node.getParent():
      return node.extra_info_index( 'etilde', self._e( node, d)), {}

    node_e = self._e( node, d)
    parent_etilde, _ = self.error( node.getParent(), d)
    return node.extra_info_index( 
            'etilde', 1/( 1/node_e + 1/parent_etilde)
           ), {} #no extra info

class Binev2013( Binev):
  def __init__( self, f, s, a=0, b=1):
    self.f = f
    self.s = s
    self.tree = t( -1, a, b)
    self.tree_hp = Tree_1D_hp( -1, a, b)
    self.error_list = []

    from error_functionals import TwoNorm
    self.errorClass = TwoNorm( self.f)

  def iteration( self):
    sorter = self.s()

    T = deepcopy( self.tree)
    for n in self.tree.leaves():
      self.__etilde( n, self.__etilde_h( n))

    b = T.branches()
    for n in b:
      p = self._p( n)
      e_p, info = self._e_p( n)
      E_hp = self._E_hp( n, T)
      print p, e_p, E_hp
      if e_p < E_hp:
        n.value( e_p)
        n.extra_info( dict( n.extra_info().items() + info.items()))
        n.trim()
        for l in self.tree.leaves():
          self.__etilde( l, self.__etilde( l) * e_p/E_hp)

    self.tree_hp = self._T_hp( T)
    sorter.add( self.tree.leaves())

    bests = sorter.find()
    return bests

  #halfway Binev 2013 page 15
  def _p( self, node):
    n = self.tree.node( *node.boundary())
    return len( n.leaves())

  def _T_hp( self, tree):
    if not tree.root().hasSameBoundary( self.tree.root()):
      raise StandardError("Illegal operation! Roots do not coincide")
      return False

    T_hp = deepcopy( tree)
    T_hp.transformTo( Tree_1D_hp)

    for l in tree.leaves():
      T_hp.node( *l.boundary()).p( self._p( l))

    return T_hp

  #TODO: improve this
  def _E_hp( self, node, T):
    leaves = self.tree.leaves()
    for node in T.leaves():
      isleaf = False
      for leaf in leaves:
        if np.allclose( node.boundary(), leaf.boundary()):
          isleaf = True
          break
      if not isleaf and node in leaves:
        leaves.remove( node)

    #leaves is now sanitized
    s = 0.0
    for n in leaves:
      e, i = self._e_p( n)
      s = s + e
    return s

  def _e_p( self, node):
    p = self._p( node)
    e, info = self.errorClass.error( node.boundary(), p)

    return e, info

  def __etilde_h( self, node):
    e1 = self._e( node, 1)

    if not node.hasParent():
      return e1

    #see Binev 2013 page 15
    return 1/(1/e1 + 1/self.__etilde_h( node.getParent()))

  def __etilde( self, node, val = False):
    return node.extra_info_index( 'etilde', val)

  def error( self, node, d = False):
    if d != False:
      raise TypeError("Argument invalid!")
      return False
    return self.__etilde( node), {}

  def needsSubdivide( self, leaf):
    return True
