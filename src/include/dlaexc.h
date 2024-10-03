      pure subroutine DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, &
     &                        WORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      logical,intent(in)     :: WANTQ
      integer,intent(in)     :: J1, LDQ, LDT, N, N1, N2
      integer,intent(out)    :: INFO
      real(wp),intent(inout) :: Q( LDQ, * ), T( LDT, * )
      real(wp),intent(out)   :: WORK( * )
      end subroutine DLAEXC
