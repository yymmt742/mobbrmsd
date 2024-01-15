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

for ai, f in zip(a, sys.argv[3:]):
    s = Path(f).stem
    #h = np.histogram2d(ai[:,0], ai[:,2], bins=(250,100), range=((0.0,2.5),(0.0,100.0)))
    h = np.histogram2d(ai[:,0], ai[:,2], bins=(250,200), range=((0.0,2.5),(0.0,10000000.0)))
    np.savetxt(s+'_hist.dat', h[0])
    plt.hist2d(ai[:,0], ai[:,2], bins=(250,200), range=((0.0,2.5),(0.0,10000000.0)))
    #plt.hist2d(ai[:,0], ai[:,2], bins=(250,64), range=((0.0,2.5),(0.0,64.0)))
    plt.savefig(s+'.png')
    plt.clf()
