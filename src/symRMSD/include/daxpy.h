pure subroutine DAXPY(N, DA, DX, INCX, DY, INCY)
  use, intrinsic :: ISO_FORTRAN_ENV, only: wp => REAL64
  real(wp), intent(in)    :: DA
  integer, intent(in)     :: INCX, INCY, N
  real(wp), intent(in)    :: DX(*)
  real(wp), intent(inout) :: DY(*)
end subroutine DAXPY

