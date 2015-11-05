from __future__ import division
from abc import ABCMeta, abstractmethod
from scipy.linalg import solve
from scipy.integrate import quad
import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt
from plotter import *

def ab2m11( a, b, x): #convert (a,b) to (-1,1) linearly
  return 2.0/(b-a) * x + (a+b)/(a-b)

class Abstract(object):
  __metaclass__ = ABCMeta

  def __init__(self, f):
    self.f = f

  """
  Assigns an error to a node; return a list [error, extra_info]
  """
  @abstractmethod
  def error(self, node, d):
    pass

class TwoNorm( Abstract):
  @staticmethod
  def norm_2( g, a, b):
    normsqr, max_error = quad( lambda x: abs(g(x)) ** 2, a, b)
    return np.sqrt(normsqr)

  @staticmethod
  def M( n, a, b):
    #based on NumAna block above (9.7), M_{jk} but from a to b.
    return np.array(
      [[(b**(k+j+1) - a**(k+j+1))/(k+j+1) for k in range(0,n)] for j in range(0,n)]
    )

  def bestpoly( self, node, d = 1):
    #d is the order which is the amount of terms in the sum (i.e. degree + 1)
    n = d
    b = np.empty((n,1))

    #see NumAna (9.7)
    M = TwoNorm.M( n, *node)
    for i in range(0, n):
      integrand = lambda x: self.f(x) * x**i
      b[i], max_error = quad(integrand, *node)

    p_n_coeffs = solve( M, b)
    p_n = lambda x: np.polyval(p_n_coeffs[::-1], x) #need to reverse order

    return p_n, p_n_coeffs

  def error( self, node, d = 1):
    p_n, p_n_coeffs = self.bestpoly( node, d)
    error = TwoNorm.norm_2( lambda x: self.f(x) - p_n(x), *node) ** 2
    #put the extra info in there
    return error, {'poly': p_n_coeffs}

class TwoNormOrth( TwoNorm):
  def __init__( self, f):
    super(self.__class__, self).__init__(f)
    from sympy.mpmath import legendre as L
    from numpy.polynomial.legendre import legval as legval
    self.L = L
    self.legval = legval

  def __gamma( self, node, n):
    a, b = node
    leg = lambda x: self.L(n, ab2m11( a, b, x)) #move
    teller, _ = quad( lambda x: self.f(x) * leg(x), a, b)
    noemer = (b-a)/(2*n+1) #see verslag
    return teller/noemer

  def bestpoly( self, node, d = 1):
    a, b = node
    leg = []
    for j in range( d):
      leg.append( self.__gamma( node, j))
    legvallambda = lambda x: self.legval( ab2m11( a, b, x), leg)
    return legvallambda, legvallambda

class Dummy(Abstract):
  def error(self, node, d = 1):
    return 2.0, None
