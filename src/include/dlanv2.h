      pure subroutine DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
      use LA_CONSTANTS, only: wp=>dp
      real(wp),intent(inout) :: A, B, C, D
      real(wp),intent(out)   :: CS, RT1I, RT1R, RT2I, RT2R, SN
      end subroutine DLANV2

