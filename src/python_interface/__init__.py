import numpy

class mobbrmsd:
  def __init__(self, d:int=3):

    if d==2:
      from .mobbrmsd_2d import driver
      self.driver = driver
    elif d==3:
      from .mobbrmsd_3d import driver
      self.driver = driver
    else:
      from .mobbrmsd import driver
      self.driver = driver
      self.driver.setup_dimension(d)

  def add_molecule(self, m=1, n=1, swp=None):

    if swp is None:
      self.driver.add_molecule(m, n, 1)
    else:
      swp_ = numpy.array(swp).reshape((-1, m))
      s = swp_.shape[0] + 1
      self.driver.add_molecule(m, n, s, swp_)

  def run(self, x:numpy.ndarray, y:numpy.ndarray, cutoff:float=float('inf'),
          difflim:float=0.0, maxeval:int=-1, rotate_y:bool=False) -> dict:

    ndim = self.driver.n_dims()
    natom = self.driver.n_atoms()

    if x.ndim == 2:
      x_ = x.transpose()
    else:
      raise IndexError

    if y.ndim == 2:
      y_ = y.transpose().reshape((-1, y.shape[1], y.shape[0]))
    elif y.ndim == 3:
      y_ = y.transpose([2,1,0])
    else:
      raise IndexError

    if x_.shape[0] != ndim or x_.shape[1] != natom or \
       y_.shape[0] != ndim or y_.shape[1] != natom:
      raise IndexError

    ntarget = y_.shape[2]
    n_header, n_int, n_float = self.driver.state_vector_lengthes()
    hret, iret, rret = self.driver.run(n_header, n_int, n_float, x_, y_, cutoff, difflim, maxeval, rotate_y)
    state = numpy.array([(i, r) for i, r in zip(iret.T, rret.T)], dtype=object)

    return {'header':hret, 'state':state}
