pure subroutine DCOPY(N, DX, INCX, DY, INCY)
  use, intrinsic :: ISO_FORTRAN_ENV, only: wp => REAL64
  integer, intent(in)   :: INCX, INCY, N
  real(wp), intent(in)  :: DX(*)
  real(wp), intent(out) :: DY(*)
end subroutine DCOPY

