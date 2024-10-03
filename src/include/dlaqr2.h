      pure subroutine DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, &
     &                        ILOZ, IHIZ, Z, LDZ, NS, ND, SR, SI, V,   &
     &                        LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK )
      use LA_CONSTANTS, only: wp=>dp
      LOGICAL,intent(in)     :: WANTT, WANTZ
      integer,intent(in)     :: IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
     &                          LDZ, LWORK, N, NH, NV, NW
      integer,intent(out)    :: ND, NS
      real(wp),intent(out)   :: SI( * ), SR( * ), T( LDT, * ), V( LDV, * ), &
     &                          WV( LDWV, * ), WORK( * )
      real(wp),intent(inout) :: H( LDH, * ), Z( LDZ, * )
      end subroutine DLAQR2
