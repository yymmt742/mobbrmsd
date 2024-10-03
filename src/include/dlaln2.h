      pure subroutine DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, &
     &                        B, LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
      use LA_CONSTANTS, only: wp=>dp
      logical,intent(in)      :: LTRANS
      integer,intent(in)      :: LDA, LDB, LDX, NA, NW
      integer,intent(out)     :: INFO
      real(wp),intent(in)     :: CA, D1, D2, SMIN, WI, WR
      real(wp),intent(out)    :: SCALE, XNORM
      real(wp),intent(out)    :: A( LDA, * ), B( LDB, * ), X( LDX, * )
      end subroutine DLALN2
