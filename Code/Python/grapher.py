import pygraphviz as gpv
import pprint
from tree import *
import os

def unique_file(file_name):
  counter = 1
  file_name_parts = os.path.splitext(file_name) # returns ('/path/file', '.ext')
  while 1:
    try:
      fd = os.open(file_name, os.O_CREAT | os.O_EXCL | os.O_RDWR)
      return file_name
    except OSError:
      pass
    file_name = file_name_parts[0] + '_' + str(counter) + file_name_parts[1]
    counter += 1

class TreeGrapher( object):
  def __init__( self, forest, outfile, debug = False, unique = False):
    self.forest = forest
    self.G = gpv.AGraph( directed=True)
    self.outfile = outfile
    self.debug = debug
    self.unique = unique

  def set_leaf( self, node):
    n = self.G.get_node( node)
    n.attr['shape'] = 'rect'

  def t( self, tree):
    if tree.hasParent():
      return tree.getParent().forest.index(tree) + 2 * self.t( tree.getParent())
    else:
      return 0

  def __treeLabel( self, tree):
    root = tree.root()
    a, b = root.boundary()
    a = round( a, 2)
    b = round( b, 2)
    t, n = [0, 2**tree.depth()] #teller, noemer

    if tree.hasParent():
      t = self.t( tree)

    strings = ["",""]
    if a != 0.:
      strings[0] = "".join( [strings[0], "%s + " % str( a)])
      strings[1] = "".join( [strings[1], "%s + " % str( a)])
    if b != 1.:
      strings[0] = "".join( [strings[0], "%s * " % str( b)])
      strings[1] = "".join( [strings[1], "%s * " % str( b)])
    string = " [ %s%i/%i, %s%i/%i ]" % ( 
      strings[0], 
      t, 
      n, 
      strings[1], 
      t+1, 
      n, 
    )
    if isinstance( tree, Tree_hp):
      string = "".join( [string, "\np=%d" % tree.p()])
    string = "".join( [string, "\ne=%s\n%s" % (tree.error(), pprint.pformat(tree.extra_info(), width=1))])

    return string

  def __recursiveGraph( self, tree, parent = None):
    self.G.add_node( self.__treeLabel( tree))
    if parent != None:
      self.G.add_edge( self.__treeLabel( parent), self.__treeLabel( tree))

    if tree.isLeaf():
      n = self.set_leaf( self.__treeLabel( tree))
    else:
      for n in tree.forest:
        self.__recursiveGraph( n, parent = tree)

  def graph( self):
    #self.G.add_edge( "boundaries\nerror value\nextra_info", "leaf")
    #self.set_leaf( "leaf")
    for tree in self.forest:
      self.__recursiveGraph( tree, parent = None)

    self.G.graph_attr.update(size="10,100")
    self.G.layout( prog = 'dot')
    if self.unique:
      self.G.draw( unique_file( self.outfile))
    else:
      self.G.draw( self.outfile)
