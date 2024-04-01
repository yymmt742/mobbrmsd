import numpy

class mobbrmsd:
  def __init__(self, d=3):
    if d==2:
      from .mobbrmsd_2d import driver
      self.driver = driver
    elif d==3:
      from .mobbrmsd_3d import driver
      self.driver = driver
    else:
      from .mobbrmsd import driver
      self.driver = driver

  def add_molecule(self, m=1, n=1, s=1, swp=None):
    if swp is None:
      self.driver.add_molecule(m, n, s)
    else:
      self.driver.add_molecule(m, n, s, numpy.array(swp))

  def run(self, x, y, cutoff=numpy.inf, difflim=0.0, maxeval=-1):
    x_ = numpy.array(x).flatten()
    y_ = numpy.array(y).flatten()
    n_ = y_.shape[0]/x_.shape[0]
    return self.driver.run(x_, y_, cutoff, difflim, maxeval, n_)
