      pure subroutine DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, &
     &                        WI, ILOZ, IHIZ, Z, LDZ, INFO )
      use LA_CONSTANTS, only: wp=>dp
      logical,intent(in)      :: WANTT, WANTZ
      integer,intent(in)      :: IHI, IHIZ, ILO, ILOZ, LDH, LDZ, N
      integer,intent(out)     :: INFO
      real(wp),intent(inout)  :: H( LDH, * ), Z( LDZ, * )
      real(wp),intent(out)    :: WI( * ), WR( * )
      end subroutine DLAHQR
