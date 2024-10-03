      pure subroutine DTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, &
     &                         LDVL, VR, LDVR, MM, M, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: HOWMNY, SIDE
      integer,intent(in)      :: LDT, LDVL, LDVR, LWORK, MM, N
      integer,intent(out)     :: M, INFO
      logical,intent(inout)   :: SELECT( * )
      real(wp),intent(inout)  :: T( LDT, * ), VL( LDVL, * ), VR( LDVR, * )
      real(wp),intent(out)    :: WORK( * )
      end subroutine DTREVC3
