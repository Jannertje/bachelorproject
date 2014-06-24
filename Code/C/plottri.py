import matplotlib.pyplot as plt
import random
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import eval_jacobi as jac
from scipy.special import eval_legendre as leg
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection, TriMesh, PathCollection, PolyCollection

def f(x, y):
  if( x + y > 1.1):
    return 1
  else: return 0

def Ginv( vol, x1, y1, x2, y2, x3, y3, x, y):
  return (1/(2.0 * vol) * ((y3-y1)*(x-x1) + (x1-x3)*(y-y1)), 1/(2.0*vol)*((y1-y2)*(x-x1) + (x2-x1)*(y-y1)));

def G( x1, y1, x2, y2, x3, y3, x, y):
  return ((x2-x1)*x + (x3-x1)*y + x1, (y2-y1)*x + (y3-y1)*y + y1)

def pkd( j, k, x, y):
  return (1-y)**j * leg(j, 2*x/(1-y) - 1) * jac( k, 2*j+1, 0, 2*y-1)

def sign(p1, p2, p3):
  return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])


def PointInAABB(pt, c1, c2):
  return c2[0] <= pt[0] <= c1[0] and \
         c2[1] <= pt[1] <= c1[1]

def PointInTriangle(pt, v1, v2, v3):
  b1 = sign(pt, v1, v2) <= 0
  b2 = sign(pt, v2, v3) <= 0
  b3 = sign(pt, v3, v1) <= 0

  return ((b1 == b2) and (b2 == b3)) and \
         PointInAABB(pt, map(max, v1, v2, v3), map(min, v1, v2, v3))

def pair_cantor( x, y):
  return (x+y)*(x+y+1)/2 + y;

draw_arrows = False
expect_coeffs = True
draw_3d = True
x = []
y = []
z = []
tris = []
coeffs = []
r = []
while True:
  line = sys.stdin.readline()
  if line:
    if line == "\n":
      break
    else:
      px, py = line.split(" ")
      x.append( float(px))
      y.append( float(py))
      z.append( 0)

while True:
  line = sys.stdin.readline()
  if line:
    if line == "\n":
      break
    else:
      p1, p2, p3 = line.split(" ")
      tris.append( [int(p1),int(p2),int(p3)])
      if expect_coeffs:
        line = sys.stdin.readline()
        stuff = line.split(" ")
        r.append( int(stuff[0]))
        coeff = []
        if stuff[0] > 0:
          for i in range( int(stuff[0])*(int(stuff[0])+1)/2):
            coeff.append( float(stuff[i+1]))
        coeffs.append( coeff)

fig = plt.figure()
ax_tri = False
if not draw_3d:
  ax_tri = fig.add_subplot(111)
else:
  ax_tri = fig.add_subplot(212)
  vols = []
  for i, t in enumerate(tris):
    vols.append( 1/2.0 * abs(-x[t[1]]*y[t[0]] + x[t[2]]*y[t[0]] + x[t[0]]*y[t[1]] - x[t[2]]*y[t[1]] - x[t[0]]*y[t[2]] + x[t[1]]*y[t[2]]));

  xs = []
  ys = []
  zs = []
  ref = []

  c=[]
  for i, t in enumerate(tris):
    if r[i] > 0:
      for _ in range( int(8000*vols[i])):
        XX = random.random()
        YY = random.random() % (1-XX)
        xf,yf = G(x[t[0]], y[t[0]], x[t[1]], y[t[1]], x[t[2]], y[t[2]], x=XX, y=YY)
        xs.append( xf)
        ys.append( yf)
        s = 0.0
        for j in range( r[i]):
          for k in range( r[i] - j):
            s += coeffs[i][pair_cantor(j, k)]*pkd( j, k, XX, YY)
        zs.append( s)#f(xf, yf))
        c.append( zs)

if draw_arrows:
  patches = []
  for t in tris:
    midpoint_edge = ((x[t[1]] + x[t[2]])/2, (y[t[1]] + y[t[2]])/2)
    arrow = mpatches.Arrow( x[t[0]], y[t[0]], -0.5 * (x[t[0]] - midpoint_edge[0]), -0.5*(y[t[0]] - midpoint_edge[1]), width=0.1)
    patches.append( arrow)
  collection = PatchCollection( patches)

  ax_tri.add_collection( collection)
ax_tri.set_aspect('equal')
#ax_tri.set_ylim([-0.2, 1.2])
plt.triplot( x, y, tris, 'g-')

if draw_3d:
  ax = fig.add_subplot(211, projection='3d')
  ax.scatter( xs, ys, zs, c=zs, cmap=plt.cm.rainbow, alpha=1.0)
  ax.set_xlim([0,1])
  ax.set_ylim([0,1])
  ax.set_zlim([0,1])
  ax.set_aspect('equal')
  """
  trian = PathCollection( TriMesh( tri.Triangulation( x, y, tris)).get_paths()).get_paths()
  polys = [ path.to_polygons()[0] for path in trian]
  print polys
  coll = PolyCollection( polys)
  ax.add_collection3d( coll, zs=0, zdir='z')
  """


plt.show()
