from .symrmsd import driver

setup_done = False

def add_molecule(m=1, n=1, s=1, swp=[0]):
  driver.set_molecule(m, n, s, swp)
  setup_done = False

def run(x, y, maxeval=None, cutoff=None):
  if maxeval != None:
    driver.maxeval = maxeval
    setup_done = False
  if cutoff != None:
    driver.cutoff = cutoff
    setup_done = False
  if setup_done == False:
    driver.setup()
  import numpy as np
  x_ = np.array(x).flatten()
  y_ = np.array(y).flatten()
  n_ = y_.shape[0]/x_.shape[0]
  return driver.run(x_, y_, n_)
