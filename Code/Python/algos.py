from abc import ABCMeta, abstractmethod

class Algo(object):
  __metaclass__ = ABCMeta

  def __init__(self, f, s, a = 0, b = 1):
    #function to approximate
    self.f = f

    #sorting class
    self.s = s

    self.a = a
    self.b = b

  """
  Return the error of this node
  """
  @abstractmethod
  def e(self, node, d):
    pass

  def iteration(self, leaves):
    sorter = self.s()

    for leaf in leaves:
      sorter.add(leaf, self.e(leaf))

    bests = sorter.find()
    return bests

  @staticmethod
  def split( leaves, bests):
    def splitNode( node):
      a, b = node
      return [ [a, (a+b)/2.0], [(a+b)/2.0, b] ]

    for (node, value) in bests:
      i = leaves.index(node)
      node1, node2 = splitNode( node)
      leaves[i] = node2
      leaves.insert( i, node1)
    
    return leaves
  
  def iterate( self):
    leaves = [ [self.a, self.b] ]
    while True:
      bests = self.iteration( leaves)
      print bests
      leaves = Algo.split( leaves, bests)


class Greedy(Algo):
  def __init__(self, f, s, d = 1):
    self.d = d
    super(Greedy, self).__init__(f, s)

    from error_functionals import TwoNorm
    self.errorClass = TwoNorm(self.f)

  def e(self, node, d = False):
    if d == False:
      d = self.d
    return self.errorClass.error(node, d)
