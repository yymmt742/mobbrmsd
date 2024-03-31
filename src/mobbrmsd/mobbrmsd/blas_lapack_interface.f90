!| Module for blas lapack interface, for general dimension.
module blas_lapack_interface
  implicit none
  private
  public  :: SGEMM, DGEMM
  public  :: D, DD, ND
  public  :: setup_dimension
!
  interface
    include 'dgemm.h'
    include 'sgemm.h'
  end interface
!
  integer, protected, save :: D = 3
  !! Spatial dimension
  integer, protected, save :: DD = 9
  !! Square spatial dimension
  integer, protected, save :: ND = 9 + 2
  !! Node size, defined by [L, G, C(D,D)]
!
contains
!
!| Sets the dimensions of the space. <br>
!  Caution, this routine affects global.
  subroutine setup_dimension(d_)
    integer, intent(in) :: d_
    D = MAX(1, d_)
    DD = D * D
    ND = DD + 2
  end subroutine setup_dimension
!
end module blas_lapack_interface

