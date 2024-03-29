!| Module for blas lapack interface, D=2.
module blas_lapack_interface
  implicit none
  private
  public  :: SGEMM, DGEMM
  public  :: D, DD, ND
  public  :: setup_dimension
!
  interface
!
    include 'dgemm.h'
    include 'sgemm.h'
!
  end interface
!
  integer, parameter :: D  = 2
  !! Spatial dimension
  integer, parameter :: DD = 4
  !! Square spatial dimension
  integer, parameter :: ND = DD + 2
  !! Node size, defined by [L, G, C(D,D)]
!
contains
!
!| Sets the dimensions of the space. <br>
!  Caution, this routine affects global.
  subroutine setup_dimension(d_)
    integer, intent(in) :: d_
  end subroutine setup_dimension
!
end module blas_lapack_interface

