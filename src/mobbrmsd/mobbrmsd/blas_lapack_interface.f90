!| Define spatial dimension, \(D\), and provide an optimized blas/lapack interface.
module blas_lapack_interface
  implicit none
  private
  public  :: D, DD, ND
  public  :: setup_dimension
  !| Spatial dimension, \(D\).
  integer, protected, save :: D = 3
  !| Square spatial dimension, \(D^2\).
  integer, protected, save :: DD = 9
  !| Node memory size, defined by \(1 + 1 + D^2\).
  !  Let \([L, G, \mathbf{C}]\) be a node,
  !  where \(L, G\in\mathbb{R}\) and \(\mathbf{C}\in\mathbb{R}^{D\times D}\).
  integer, protected, save :: ND = 9 + 2
contains
!| Sets the dimensions of the space. <br>
!  Caution, this routine affects global.
  subroutine setup_dimension(d_)
    integer, intent(in) :: d_
    D = MAX(1, d_)
    DD = D * D
    ND = DD + 2
  end subroutine setup_dimension
end module blas_lapack_interface

