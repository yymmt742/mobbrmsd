import numpy as np
import sys
from pathlib import Path
import matplotlib.pyplot as plt

for f in sys.argv[1:]:
    print(f)
    a = np.loadtxt(f)
    s = Path(f).stem
    h = np.histogram2d(a[:,0], a[:,2], bins=(250,100), range=((0.0,2.5),(0.0,100.0)))
    np.savetxt(s+'_hist.dat', h[0])
    plt.hist2d(a[:,0], a[:,2], bins=(250,64), range=((0.0,2.5),(0.0,64.0)))
    plt.savefig(s+'.png')
    plt.clf()
