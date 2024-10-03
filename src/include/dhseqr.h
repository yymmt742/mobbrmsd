      pure subroutine DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, &
     &                        Z, LDZ, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: COMPZ, JOB
      integer,intent(in)      :: IHI, ILO, LDH, LDZ, LWORK, N
      integer,intent(out)     :: INFO
      real(wp),intent(inout)  :: H( LDH, * ), Z( LDZ, * )
      real(wp),intent(out)    :: WI( * ), WORK( * ), WR( * )
      end subroutine DHSEQR
