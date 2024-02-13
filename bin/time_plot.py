import sys
import numpy as np
import matplotlib.pyplot as plt

c = ["#FF3300", "#FF9900", "#339966", "#336699"]
t = ['0.1', '0.2', '0.3', '0.4']
n = np.loadtxt(sys.argv[1])[:,0]
a = np.array( [np.loadtxt(p)[:,2] for p in sys.argv[1:]])
la = (np.mean(np.log10(a[:,-1]))-np.mean(np.log10(a[:,5])))/(n[-1]-n[5])
lb = np.mean(np.log10(a[:,-1])) - la * 16
print('slope =', la, 10**la)
exit()
ax = plt.axes()
ax.set_aspect(0.3535*16/7)
ax.set_ylim(-5,2)
ax.plot(n, n*la+lb-1.0, label='slope', color='black', linestyle='--')
for ai, ti, ci in zip(a, t, c):
  ax.plot(n, np.log10(ai), label=ti, color=ci)
ax.legend()
plt.savefig('timing.png')
plt.savefig('timing.eps')
