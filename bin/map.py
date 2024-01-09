import numpy as np
import sys
from pathlib import Path
import matplotlib.pyplot as plt

s = Path(sys.argv[1]).stem

a = np.loadtxt(sys.argv[1])
for f in sys.argv[2:]:
    a += np.loadtxt(f)

a /= len(sys.argv[1:])

plt.imshow(a)
plt.savefig(s+'.png')
plt.clf()
