pure subroutine DCOPY(N, X, INCX, Y, INCY)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL64
  integer, intent(in)   :: INCX, INCY, N
  real(rk), intent(in)  :: X(*)
  real(rk), intent(out) :: Y(*)
end subroutine DCOPY

