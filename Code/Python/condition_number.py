import math
import numpy as np
from numpy.linalg import cond
import matplotlib.pyplot as plt
import sys

def integral( xpow, ypow):
  return (math.factorial( xpow) * math.factorial( ypow) + 0.0)/(math.factorial( 2 + xpow + ypow))

def pair( x, y):
  return (x+y)*(x+y+1)/2 + y

def unpair( z):
  w = int( math.floor( (math.sqrt(8*z + 1) - 1)/2))
  t = (w*w + w)/2
  return (w - z+t, z - t)

N = 10
X = []
for n in range(N):
  numbas = (n+2)*(n+1)/2
  H = np.zeros((numbas,numbas))

  for i in range(0, numbas):
    for j in range(0, numbas):
      a, b = unpair( i)
      c, d = unpair( j)
      H[i, j] = integral( a + c, b + d)

  X.append(cond(H))

print X
Y = [((n+2)*(n+1) +1)*((n+2)*(n+1)/2+1) for n in range(N)]
Z = [(n+2)*(n+1)/2 for n in range(N)]
print Y
plt.plot( range(N), X)
plt.plot( range(N), Y)
plt.legend(["Monomial basis", "PKD basis"])
plt.yscale('log')
plt.xlabel('Highest degree of basis elements')
plt.ylabel('Condition number of basis')
plt.show()
