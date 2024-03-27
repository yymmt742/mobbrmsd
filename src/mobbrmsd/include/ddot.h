pure function DDOT(N, X, INCX, Y, INCY)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL64
  integer, intent(in)  :: INCX, INCY, N
  real(rk), intent(in) :: X(*), Y(*)
  real(rk)             :: DDOT
end function DDOT

