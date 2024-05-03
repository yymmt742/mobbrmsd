import numpy as np
import sys
from pathlib import Path
import matplotlib.pyplot as plt

n = int(sys.argv[1])
m = int(sys.argv[2])
#n = int(np.sqrt(len(sys.argv[1:])))
#m = int(np.sqrt(len(sys.argv[1:])))
a = np.array([np.loadtxt(f) for f in sys.argv[3:]])
mean = np.mean(a[:,:,2], 1).reshape([n,m])
var = np.var(a[:,:,2], 1).reshape([n,m])

np.savetxt('mean.dat', mean)
np.savetxt('var.dat', var)

hmin = 0
hmax = 600000

for ai, f in zip(a, sys.argv[3:]):
    s = Path(f).stem
    h = np.histogram2d(ai[:,0], ai[:,2], bins=(250,200), range=((0.0,2.5),(hmin,hmax)))
    np.savetxt(s+'_hist.dat', h[0])
    plt.hist2d(ai[:,0], ai[:,2], bins=(250,200), range=((0.0,2.5),(hmin,hmax)))
    plt.savefig(s+'.png')
    plt.clf()
