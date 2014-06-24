import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
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

fig = plt.figure()
plt.triplot( x, y, tris, 'go-')
plt.show()
