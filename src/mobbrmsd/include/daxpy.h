pure subroutine DAXPY(N, A, X, INCX, Y, INCY)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL64
  real(rk), intent(in)    :: A
  integer, intent(in)     :: INCX, INCY, N
  real(rk), intent(in)    :: X(*)
  real(rk), intent(inout) :: Y(*)
end subroutine DAXPY

