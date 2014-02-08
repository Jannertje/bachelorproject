from abc import ABCMeta, abstractmethod

class Algo(object):
  __metaclass__ = ABCMeta

  def __init__(self, f, s):
    #function to approximate
    self.f = f

    #sorting class
    self.s = s

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
    print bests

class Greedy(Algo):
  def __init__(self, f, s):
    super(Greedy, self).__init__(f, s)

    from error_functionals import Dummy
    self.errorClass = Dummy(self.f)

  def e(self, node, d = 1):
    return self.errorClass.error(node, d)
