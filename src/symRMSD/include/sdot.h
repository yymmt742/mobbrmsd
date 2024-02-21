pure function SDOT(N, X, INCX, Y, INCY)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL32
  integer, intent(in)  :: INCX, INCY, N
  real(rk), intent(in) :: X(*), Y(*)
  real(rk)             :: SDOT
end function SDOT

