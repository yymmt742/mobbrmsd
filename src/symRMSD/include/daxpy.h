pure subroutine DAXPY(N,DA,DX,INCX,DY,INCY)
  use mod_params, only :  wp=>R8
  real(wp), intent(in)    :: DA
  integer,  intent(in)    :: INCX,INCY,N
  real(wp), intent(inout) :: DX(*),DY(*)
end subroutine DAXPY
