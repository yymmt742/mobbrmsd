pure function DDOT(N, DX, INCX, DY, INCY)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL64
  integer, intent(in)      :: N
  real(rk), intent(in)     :: DX(*)
  integer, intent(in)      :: INCX
  real(rk), intent(in)     :: DY(*)
  integer, intent(in)      :: INCY
  real(rk)                 :: DDOT
end function DDOT

