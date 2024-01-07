import numpy as np
import sys
import matplotlib.pyplot as plt

n = int(np.sqrt(len(sys.argv[1:])))
a = np.array([np.loadtxt(f) for f in sys.argv[1:]])
mean = np.mean(a[:,:,2], 1).reshape([n,n])
var = np.var(a[:,:,2], 1).reshape([n,n])

np.savetxt('mean.dat', mean)
np.savetxt('var.dat', var)

plt.imshow(mean)
plt.show()
