      pure subroutine DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, &
     &                        SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, &
     &                        LDU, NV, WV, LDWV, NH, WH, LDWH )
      use LA_CONSTANTS, only: wp=>dp
      LOGICAL,intent(in)     :: WANTT, WANTZ
      integer,intent(in)     :: IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, &
     &                          LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
      real(wp),intent(inout) :: H( LDH, * ), Z( LDZ, * ), SR( * ), SI( * )
      real(wp),intent(out)   :: U( LDU, * ), V( LDV, * ), WH( LDWH, * ), WV( LDWV, * )
      end subroutine DLAQR5
