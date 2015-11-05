import matplotlib.pyplot as plt
from scipy.special import eval_legendre as leg
import numpy as np
import sys
from math import ceil

def f( x):
  if( x < 3): return np.exp(x)/np.exp(3);
  else: return np.sin(np.pi *x/2);

def pol( a, b, x, coeff):
  s = 0
  for i in range( len( coeff)):
    s = s + coeff[i] * leg( i, 2.0/(b-a) * x + (a+b)/(a-b))
  return s

coeffs = []
r = []
a = []
b = []

while True:
  line = sys.stdin.readline()
  if line:
    if line == "\n":
      break
    else:
      array = line.split(" ")
      cura = float(array[0])
      curb = float(array[1])
      curr = int(array[2])
      a.append( cura)
      b.append( curb)
      r.append( curr)
      curcoeffs = []
      for i in range( curr):
        curcoeffs.append( float( array[i+3]))
      coeffs.append( curcoeffs)

plt.subplot(121)
XX = np.linspace( min(a), max(b), 10000)
plt.plot( XX, [f(x) for x in XX])

print a, b, r, coeffs
maxs = []
mins = []
for i in range( len(coeffs)):
  X = np.linspace( a[i], b[i], 10000*(b[i]-a[i]))
  if len(X):
    pX = np.array([pol(a[i], b[i], x, coeffs[i]) for x in X])
    maxs.append( max( pX))
    mins.append( min( pX))
    plt.plot( X, pX)
#plt.ylim( min(mins) - 0.2, max(maxs) + 0.2)

plt.subplot(122)
for i in range( len(coeffs)):
  plt.bar( a[i], [r[i]], b[i]-a[i], color=str(1 - 1.0/(1+r[i])))

plt.show()
