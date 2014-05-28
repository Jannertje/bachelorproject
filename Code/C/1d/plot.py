import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
f = open( "test.txt", "r")
nbins = 100000
bins = [0]*nbins
minbins = [2.0]*nbins
maxbins = [-2.0]*nbins
for line in f:
  if line[0] == "I":
    continue
  x, y = line.split(" ")
  x = float(x)
  y = float(y)
  curbin = int((x+1)*nbins/2.0)
  bins[curbin] += 1
  if minbins[curbin] > y:
    minbins[curbin] = y
  if maxbins[curbin] < y:
    maxbins[curbin] = y
f.close()

for i in range( nbins):
  if( bins[i] == 0):
    maxbins[i] = 0
    minbins[i] = 0

plt.subplot(211)
plt.plot( bins)
plt.subplot(212)
plt.plot( np.array(maxbins) - np.array(minbins))
plt.show()
