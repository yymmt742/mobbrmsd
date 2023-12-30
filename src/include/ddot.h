pure function DDOT(N,DX,INCX,DY,INCY)
  use mod_params, only :  wp=>R8
  integer,intent(in)  :: INCX,INCY,N
  real(wp),intent(in) :: DX(*),DY(*)
  real(wp)            :: DDOT
end function DDOT
