      pure function DDOT(N,DX,INCX,DY,INCY)
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)  :: INCX,INCY,N
      real(wp),intent(in) :: DX(*),DY(*)
      real(wp)            :: DDOT
      end function DDOT
