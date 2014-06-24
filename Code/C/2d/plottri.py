import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

get_input = True
draw_arrows = True
if get_input:
  leaves = []
  while True:
    line = sys.stdin.readline()
    if line:
      if line == "\n":
        break
      else:
        ax, ay, bx, by, cx, cy = line.split(" ")
        leaves.append( ((float(ax), float(ay)), (float(bx), float(by)), (float(cx), float(cy))))

#leaves = [ [(-1.0, -1.0), (1, -1), (-1, 1)]]
"""
#tegenvb
s3 = np.sqrt(3.0)
leaves = [ [(-1.0, -1/s3), (0.0, 0.0), (0.0, 2/s3)],
           [(1.0, -1/s3), (0.0, 0.0), (-1.0, -1/s3)],
           [(0.0, 2/s3), (0.0, 0.0), (1.0, -1/s3)]]
           """
#leaves = [ [(0.0, 0.0), (1.0,0.0), (0.0,1.0)]]
leaves = np.array(leaves)
points = list(set( [ tuple(i) for i in leaves.reshape((1,-1,2))[0].tolist()]))
tris = []
for leaf in leaves:
  tris.append( [points.index(tuple(leaf[i])) for i in range(3)])
x = [p[0] for p in points]
y = [p[1] for p in points]
#t = tri.Triangulation( [p[0] for p in points], [p[1] for p in points], tris)

fig, ax = plt.subplots()
plt.axis('equal')
#plt.axis('off')
if draw_arrows:
  patches = []
  for l in leaves:
    midpoint_edge = ((l[1][0] + l[2][0])/2, (l[1][1] + l[2][1])/2)
    print l
    print midpoint_edge
    arrow = mpatches.Arrow( l[0][0], l[0][1], -0.5 * (l[0][0] - midpoint_edge[0]), -0.5*(l[0][1] - midpoint_edge[1]), width=0.1)
    patches.append( arrow)
  collection = PatchCollection( patches)

  ax.add_collection( collection)
plt.triplot( x, y, tris, 'go-')
plt.show()
