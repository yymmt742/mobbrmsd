
class mobbrmsd:
  def __init__(self, d=3):
    import numpy
    self.to_numpy = numpy.array
    if d==2:
      from .mobbrmsd_2d import driver
      self.driver = driver
    elif d==3:
      from .mobbrmsd_3d import driver
      self.driver = driver
    else:
      from .mobbrmsd import driver
      self.driver = driver

  def add_molecule(self, m=1, n=1, swp=None):
    if swp is None:
      self.driver.add_molecule(m, n, 1)
    else:
      swp_ = self.to_numpy(swp).reshape((-1, m))
      s = swp_.shape[0] + 1
      self.driver.add_molecule(m, n, s, swp_)

  def run(self, x, y, cutoff=float('inf'), difflim=0.0, maxeval=-1):
    x_ = self.to_numpy(x).flatten()
    y_ = self.to_numpy(y).flatten()
    n_ = y_.shape[0]/x_.shape[0]
    return self.driver.run(x_, y_, cutoff, difflim, maxeval, n_)
