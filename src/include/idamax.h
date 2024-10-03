      pure function IDAMAX(N,DX,INCX)
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)  :: INCX,N
      real(wp),intent(in) :: DX(*)
      integer             :: IDAMAX
      end function IDAMAX
