import numpy as np

def poly( p, x, y):
  Nx = len(x)
  Np = (p+1)*(p+2)/2
  s = -np.ones(Nx)
  t = y
  dex = np.abs(y-1) > 1e-10
  s[dex] = 2*(1+x[dex])/(1-y[dex])-1

  V = np.zeros((Nx,Np))
  ll = 0
  tfact = np.ones(Nx)

poly(1, [0], [0])
