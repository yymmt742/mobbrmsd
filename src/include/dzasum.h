      pure function DZASUM(N,ZX,INCX)
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     :: INCX,N
      real(wp)               :: DZASUM
      complex(wp),intent(in) :: ZX(*)
      end function DZASUM
