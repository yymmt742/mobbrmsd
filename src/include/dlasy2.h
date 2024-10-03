      pure subroutine DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, &
     &                        TR, LDTR, B, LDB, SCALE, X, LDX, XNORM, &
     &                        INFO )
      use LA_CONSTANTS, only: wp=>dp
      logical,intent(in)   :: LTRANL, LTRANR
      integer,intent(in)   :: ISGN, LDB, LDTL, LDTR, LDX, N1, N2
      integer,intent(out)  :: INFO
      real(wp),intent(out) :: SCALE, XNORM
      real(wp),intent(in)  :: B( LDB, * ), TL( LDTL, * ), TR( LDTR, * )
      real(wp),intent(out) :: X( LDX, * )
      end subroutine DLASY2
