pure function SDOT(N, DX, INCX, DY, INCY)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL32
  integer, intent(in)      :: N
  real(rk), intent(in)     :: DX(*)
  integer, intent(in)      :: INCX
  real(rk), intent(in)     :: DY(*)
  integer, intent(in)      :: INCY
  real(rk)                 :: SDOT
end function SDOT

