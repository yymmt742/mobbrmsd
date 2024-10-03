      pure function IZAMAX(N,ZX,INCX)
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     :: INCX,N
      COMPLEX(wp),intent(in) :: ZX(*)
      integer                :: IZAMAX
      end function IZAMAX
