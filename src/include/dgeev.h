      pure subroutine DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, &
     &                       LDVL, VR, LDVR, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: JOBVL, JOBVR
      integer,intent(in)      :: LDA, LDVL, LDVR, LWORK, N
      integer,intent(out)     :: INFO
      real(wp),intent(inout)  :: A( LDA, * )
      real(wp),intent(out)    :: VL( LDVL, * ), VR( LDVR, * ), &
     &                           WI( * ), WORK( * ), WR( * )
      end subroutine DGEEV
