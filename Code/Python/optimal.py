from error_functionals import TwoNorm
from copy import deepcopy
from tree import Tree_1D as t
from grapher import TreeGrapher

class OptimalTree( object):
  def __init__( self, f, a = 0, b = 1, d = 1, graph = 0):
    self.f = f
    self.trees = [t(-1, a, b)]
    self.d = d
    self.graph = graph
    self.errorClass = TwoNorm( self.f)

  def iteration( self):
    newtrees = []
    for tree in self.trees:
      for leaf in tree.leaves():
        newtree = deepcopy(tree)
        newtree.subdivide( leaf)
        newtrees.append( newtree)

    self.trees = newtrees
    optimal = self.trees[0]
    optimalError = self.totalError( optimal)
    for tree in self.trees:
      tmpError = self.totalError( tree)
      if tmpError < optimalError:
        optimal = tree
        optimalError = tmpError

    print optimalError
    return (optimal, optimalError)

  def iterate( self, n):
    optimalTrees = []
    for i in range( n):
      optimalTrees.append( self.iteration())
      if self.graph:
        TreeGrapher( [optimalTrees[-1][0]], 'optimalTree.pdf', unique = False).graph()

    return optimalTrees

  def totalError( self, tree):
    som = 0

    for leaf in tree.leaves():
       e, info = self.errorClass.error( leaf.boundary(), self.d)
       leaf.value( e)
       leaf.extra_info( dict( leaf.extra_info().items() + info.items()))
       som = som + e

    return som
