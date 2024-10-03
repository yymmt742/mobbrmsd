      pure function DZNRM2( n, x, incx )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)        :: incx, n
      complex(wp),intent(in)    :: x(*)
      real(wp)                  :: DZNRM2
      end function DZNRM2
