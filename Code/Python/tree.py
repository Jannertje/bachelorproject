from abc import ABCMeta, abstractmethod
import numpy as np
from collections import deque

class Tree( object):
  __metaclass__ = ABCMeta
  K = 0

  def __init__( self, value, a, b, forest = None, parent = None, extra_info = None):
    self.__value = value
    self.__error = -1
    self.__extra_info = {}
    self.__parent = parent
    self.a = a
    self.b = b
    self.forest = []

    if isinstance( forest, list):
      if len(forest) > K:
        raise IndexError("Too many children!")
      for tree in forest:
        if isinstance( tree, self.__class__):
          self.forest.append( tree.setParent( self))
        else:
          try:
            val, a, b = tree
            if type( val) == float:
              self.forest.append( self.__class__( val, a, b, parent=self))
          except TypeError:
            raise TypeError("Type not supported!")

  def hasSameBoundary( self, node):
    return np.allclose( self.boundary(), node.boundary())

  def transformTo( self, class_):
    assert( issubclass( class_, Tree))
    self.__class__ = class_
    for n in self.forest:
      n.transformTo( class_)

  def isLeaf( self, tree = None):
    if tree == None:
      tree = self
    return len(tree.forest) == 0

  def height( self, tree = None):
    if tree == None:
      tree = self

    if tree.isLeaf():
      return 0
    else:
      return 1 + max( [n.height() for n in tree.forest])

  def trim( self):
    self.forest = []

  #return the branches in a fine-to-coarse manner
  #reversed breadth-first
  #TODO: this returns the whole tree. see C implementation.
  def branches( self):
    lst = []
    q = deque([])
    q.append( self)
    while q:
      n = q.popleft()
      lst.append( n)
      for l in n.forest:
        q.append( l)

    return reversed( lst)

  #finds the node with boundary (a,b) as
  def node( self, a, b):
    a2, b2 = self.boundary()
    if np.allclose( [a, b], [a2, b2]):
      return self

    for n in self.forest:
      ret = n.node( a, b)
      if ret:
        return ret

    return False

  def depth( self, tree = None):
    if tree == None:
      tree = self

    if tree.hasParent():
      return 1 + tree.getParent().depth()
    else:
      return 0

  def root( self):
    if self.hasParent():
      return self.getParent().root()
    else:
      return self

  def setParent( self, par):
    if par:
      self.__parent = par

    return self

  def hasParent( self):
    return self.getParent() != None

  def getParent( self):
    return self.__parent

  def siblings( self):
    if not self.getParent():
      raise Error("No parent!")

    return self.getParent().forest

  def sum_of_leaves( self):
    s = [node.error() for node in self.leaves()]
    return sum( s)

  def error( self, val = False):
    if val != False:
      self.__error = val

    return self.__error

  def value( self, val = False):
    if val != False:
      self.__value = val

    return self.__value

  def extra_info_index( self, index, val = None):
    if not self.__extra_info:
      self.__extra_info = {}

    if val != None:
      self.__extra_info[index] = val
      return self.__extra_info[index]

    if index in self.__extra_info:
      return self.__extra_info[index]
    else:
      raise IndexError("index not set!")
      print "WTF BRUH"
      return False

  def extra_info( self, info = None):
    if info != None:
      self.__extra_info = info

    return self.__extra_info

  def boundary( self):
    return [self.a, self.b]

  def leaves( self):
    res = []
    for node in self.forest:
      res.extend( node.leaves())

    if len(res) == 0:
      res = [self]

    return res


  @abstractmethod
  def subdivide( self, node):
    pass
      

class Tree_1D( Tree):
  K = 2

  def subdivide( self, node = False):
    if node == False:
      node = self
    else:
      node = self.node( node.a, node.b)
      if node == False:
        raise IndexError("node not found in current tree!")
        return False

    if len(node.forest) > 0:
      print "subdividing non-leaf...beware"
    
    a = node.a
    b = node.b
    h = (a+b)/2.0

    node.forest = [ node.__class__( -1., a, h, parent = node),
                    node.__class__( -1., h, b, parent = node) ]

    return node.forest

class Tree_hp( Tree):
  _p = -1

  def p( self, val = False):
    if val != False:
      print "p set to", val
      self._p = val

    return self._p

class Tree_1D_hp( Tree_hp, Tree_1D):
  pass
