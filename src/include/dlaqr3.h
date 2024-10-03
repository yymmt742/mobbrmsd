      pure subroutine DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, &
     &                        ILOZ, IHIZ, Z, LDZ, NS, ND, SR, SI, V,  &
     &                        LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK )
      use LA_CONSTANTS, only: wp=>dp
      logical,intent(in)     :: WANTT, WANTZ
      integer,intent(in)     :: IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
     &                          LDZ, LWORK, N, NH, NV, NW
      integer,intent(out)    :: ND, NS
      real(wp),intent(inout) :: H( LDH, * ), Z( LDZ, * )
      real(wp),intent(out)   :: SI( * ), SR( * ), V( LDV, * ), T( LDT, * ), &
     &                          WORK( * ), WV( LDWV, * )
      end subroutine DLAQR3
