import pygraphviz as gpv
import pprint

class TreeGrapher( object):
  def __init__( self, forest, outfile, debug = True):
    self.forest = forest
    self.G = gpv.AGraph( directed=True)
    self.outfile = outfile
    self.debug = debug

  def set_leaf( self, node):
    n = self.G.get_node( node)
    n.attr['shape'] = 'rect'

  def t( self, tree):
    if tree.hasParent():
      return tree.getParent().forest.index(tree)
    else:
      return 0

  def __treeLabel( self, tree):
    root = tree.root()
    a, b = root.boundary()
    a = round( a, 2)
    b = round( b, 2)
    t, n = [0, 2**tree.depth()] #teller, noemer

    if tree.hasParent():
      t = self.t( tree) + 2 * self.t( tree.getParent())

    strings = ["",""]
    if a != 0.:
      strings[0] = "".join( [strings[0], "%s + " % str( a)])
      strings[1] = "".join( [strings[1], "%s + " % str( a)])
    if b != 1.:
      strings[0] = "".join( [strings[0], "%s * " % str( b)])
      strings[1] = "".join( [strings[1], "%s * " % str( b)])
    string = " [ %s%i/%i, %s%i/%i ] \n %s" % ( 
      strings[0], 
      t, 
      n, 
      strings[1], 
      t+1, 
      n, 
      tree.value()
    )
    if self.debug:
      string = "".join( [string, "\n%s" % pprint.pformat(tree.extra_info(), width=1)])

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
    self.G.add_edge( "boundaries\nerror value\nextra_info", "leaf")
    self.set_leaf( "leaf")
    for tree in self.forest:
      self.__recursiveGraph( tree, parent = None)

    self.G.graph_attr.update(size="10,100")
    self.G.layout( prog = 'dot')
    self.G.draw( self.outfile)
