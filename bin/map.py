import numpy as np
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

s = Path(sys.argv[1]).stem

a = np.loadtxt(sys.argv[1])
#for f in sys.argv[2:]:
  #a += np.loadtxt(f).T

fig, ax = plt.subplots()

amin = 0.2
amax = 0.8
#hmin, hmax = 25, 46
hmin, hmax = 300, 600
hmap = ax.imshow(a, extent=(amin,amax,1.0,0.0), cmap='cividis')
ax.set_xticks([0.2,0.4,0.6,0.8], minor=False)
ax.set_yticks(np.linspace(0.0, 1.0, 6), minor=False)
ax.xaxis.tick_top()
hmap.set_clim(hmin,hmax)
cbar = plt.colorbar(hmap, ax=ax)

plt.savefig(s+'.png')
plt.savefig(s+'.eps')
plt.clf()
