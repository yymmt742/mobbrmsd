!| Module for blas lapack interface, D=3.
module blas_lapack_interface
  implicit none
  private
  public  :: SCOPY, DCOPY
  public  :: SGEMM, DGEMM
  public  :: D, DD, ND
  public  :: setup_dimension
!
  interface
!
    include 'dcopy.h'
    include 'scopy.h'
!
    include 'dgemm.h'
    include 'sgemm.h'
!
  end interface
!
  integer, parameter :: D  = 3
  !! Spatial dimension
  integer, parameter :: DD = 9
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

