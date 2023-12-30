pure subroutine DCOPY(N,DX,INCX,DY,INCY)
  use mod_params, only :  wp=>R8
  integer,intent(in)   :: INCX,INCY,N
  real(wp),intent(in)  :: DX(*)
  real(wp),intent(out) :: DY(*)
end subroutine DCOPY
