import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np
from random import random
import sys

leaves = np.array(
    [[(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)], 
     [(0.0, 1.0), (1.0, 0.0), (0.9, 0.9)]]
    )
points = list(set( [ tuple(i) for i in leaves.reshape((1,-1,2))[0].tolist()]))
tris = []
for leaf in leaves:
  tris.append( [points.index(tuple(leaf[i])) for i in range(3)])
x = [p[0] for p in points]
y = [p[1] for p in points]
print tris
#t = tri.Triangulation( [p[0] for p in points], [p[1] for p in points], tris)

def f(x,y):
  return np.sin(x)*np.exp(y)

X = np.linspace(0,10, 25)
Y = np.linspace(0,5, 25)
X, Y = np.meshgrid( X, Y)
Z = f(X,Y)

print X, Y, Z

fig = plt.figure()
ax = Axes3D( fig)
ax.plot_wireframe( X, Y, Z)
plt.show()
plt.triplot( x, y, tris, 'go-')
plt.show()
