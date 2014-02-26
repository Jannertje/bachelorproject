from __future__ import division
from abc import ABCMeta, abstractmethod
from scipy.linalg import solve
from scipy.integrate import quad
import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt
from plotter import *

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
  def __init__( self, f):
    super(TwoNorm, self).__init__( f)

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

  def error( self, node, d = 1):
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

    error = TwoNorm.norm_2( lambda x: self.f(x) - p_n(x), *node) ** 2
    #put the extra info in there
    return error, {'poly': p_n}

class Dummy(Abstract):
  def error(self, node, d = 1):
    return 2.0, None
