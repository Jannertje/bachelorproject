from abc import ABCMeta, abstractmethod
class Abstract(object):
  __metaclass__ = ABCMeta

  """
  Adds a (node,value) pair to the sorting algorithm
  """
  @abstractmethod
  def add(self, node, value):
    pass

  """
  Find the (possibly many) highest values in the list
  """
  @abstractmethod
  def find(self):
    pass

class Histogram(Abstract):
  pass

class Sort(Abstract):
  def __init__(self):
    self.values = []

  def add(self, node, value):
    self.values.append((node,value))

  def find(self):
    bests = []
    highest = -1.0
    
    for (node,value) in self.values:
      if value > highest:
        bests = [(node,value)]
        highest = value
      elif value == highest:
        bests.append((node,value))

    return bests
