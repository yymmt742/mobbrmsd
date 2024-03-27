from .mobbrmsd import driver
import numpy as np

def add_molecule(m=1, n=1, s=1, swp=None):
  if swp is None:
    driver.add_molecule(m, n, s)
  else:
    swp_ = np.array(swp)
    driver.add_molecule(m, n, s, swp_)

def run(x, y, cutoff=np.inf, difflim=0.0, maxeval=-1):
  x_ = np.array(x).flatten()
  y_ = np.array(y).flatten()
  n_ = y_.shape[0]/x_.shape[0]
  return driver.run(x_, y_, cutoff, difflim, maxeval, n_)
