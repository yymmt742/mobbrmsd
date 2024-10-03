      pure subroutine DRSCL( N, SA, SX, INCX )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)      :: INCX, N
      real(wp),intent(in)     :: SA
      real(wp),intent(inout)  :: SX( * )
      end subroutine DRSCL
