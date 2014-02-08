from abc import ABCMeta, abstractmethod

class Abstract(object):
  __metaclass__ = ABCMeta

  def __init__(self, f):
    self.f = f

  """
  Assigns an error to a node
  """
  @abstractmethod
  def error(self, node, d):
    pass

class Dummy(Abstract):
  def error(self, node, d = 1):
    return 2.0
