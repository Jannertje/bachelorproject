from abc import ABCMeta, abstractmethod
import numpy as np
import math

class Abstract(object):
  __metaclass__ = ABCMeta

  """
  Adds a (node,value) pair to the sorting algorithm
  """
  @abstractmethod
  def add(self, leaf):
    pass

  """
  Find the (possibly many) highest values in the list
  """
  @abstractmethod
  def find(self):
    pass

class Histogram(Abstract):
  def __init__( self):
    self.bins = {}

  def add( self, leaf):
    if isinstance( leaf, list):
      for l in leaf:
        self.add( l)
    else:
      def power_two(n):
        return int( math.log( n, 2))
      
      p = power_two( leaf.value())
      if p in self.bins:
        self.bins[p].append( leaf)
      else:
        self.bins[p] = [leaf]

  def find( self):
    m = max( self.bins.keys(), key=int)
    return self.bins[m]

class Sort(Abstract):
  def __init__(self):
    self.values = []

  def add(self, leaf):
    if isinstance( leaf, list):
      for l in leaf:
        self.add( l)
    else:
      self.values.append( leaf)

  def find(self):
    bests = []
    highest = -1.0
    
    for leaf in self.values:
      value = leaf.value()
      if value > highest:
        bests = [leaf]
        highest = value
      elif value == highest:
        bests.append( leaf)

    return bests
