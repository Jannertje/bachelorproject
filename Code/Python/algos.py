from abc import ABCMeta, abstractmethod
from tree import Tree_1D as t
from tree import Tree_1D_hp as t_hp
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

class Binev2004( Algo):
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

class Binev2004First( Binev2004):
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

class Binev2004Second( Binev2004):
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

class Binev2007( Binev2004):
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

class Binev2013( Algo):
  def __init__( self, f, s, a=0, b=1):
    self.f = f
    self.s = s
    self.tree = t( -1, a, b)
    self.tree_hp = t_hp( {}, a, b)

    from error_functionals import TwoNorm
    self.errorClass = TwoNorm( self.f)

  def interation( self):
    sorter = self.s()

  #halfway Binev 2013 page 15
  def _p( self, node):
    return len( node.leaves())

  def _T_hp( self, tree):
    if tree.root() != self.tree.root():
      raise Error("Illegal operation! Roots do not coincide")
      return False

    T_hp = t_hp( {}, *self.tree.boundary())
    T_hp.copy_from( self.tree)
    if tree.isLeaf():
      T_hp.node( *tree.boundary()).p( 1)
      return T_hp

    for l in tree.leaves():
      T_hp.node( *l.boundary()).p( self._p( l))
      return T_hp

  def _E_hp( self, node, tree):
    pass

  def _e_p( self, node, p):
    v = node.value()
    if p in v:
      return v[p]

    e, info = self.errorClass.error( node.boundary(), p)
    node.value( e, p)
    node.extra_info( dict( node.extra_info().items() + info.items()))
    return e

  def __etilde_h( self, node):
    e1 = self._e( node, 1)

    if not node.hasParent():
      return e1

    #see Binev 2013 page 15
    return 1/(1/e1 + 1/self.__etilde_h( node.getParent()))

  def error( self, node, d):
    return 0.
