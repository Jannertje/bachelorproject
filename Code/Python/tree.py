from abc import ABCMeta, abstractmethod

class Tree( object):
  __metaclass__ = ABCMeta
  K = 0

  def __init__( self, value, a, b, forest = None, parent = None, extra_info = None):
    self.__value = value
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

  def setParent( self, par):
    if par:
      self.__parent = par

    return self

  def getParent( self):
    return self.__parent

  def siblings( self):
    if not self.getParent():
      raise Error("No parent!")

    return self.getParent().forest

  def sum_of_leaves( self):
    return sum( [node.value() for node in self.leaves()])

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
      return False

  def extra_info( self, info = None):
    if info != None:
      self.__extra_info = info

    return self.__extra_info

  def boundary( self):
    return [self.a, self.b]

  @abstractmethod
  def leaves( self):
    pass

  @abstractmethod
  def subdivide( self):
    pass
      

class Tree_1D( Tree):
  K = 2

  def leaves( self):
    res = []
    for node in self.forest:
      res.extend( node.leaves())

    if len(res) == 0:
      res = [self]

    return res

  def subdivide( self):
    if len(self.forest) > 0:
      print "subdividing non-leaf...beware"
    
    a = self.a
    b = self.b
    h = (a+b)/2.0

    self.forest = [ self.__class__( -1., a, h, parent = self),
                    self.__class__( -1., h, b, parent = self) ]

    return self.forest
