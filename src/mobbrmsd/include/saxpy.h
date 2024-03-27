pure subroutine SAXPY(N, A, X, INCX, Y, INCY)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL32
  real(rk), intent(in)    :: A
  integer, intent(in)     :: INCX, INCY, N
  real(rk), intent(in)    :: X(*)
  real(rk), intent(inout) :: Y(*)
end subroutine SAXPY

