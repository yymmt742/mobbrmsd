      pure subroutine DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, &
     &                        WORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: COMPQ
      integer,intent(in)      :: IFST, ILST, LDQ, LDT, N
      integer,intent(out)     :: INFO
      real(wp),intent(inout)  :: Q( LDQ, * ), T( LDT, * )
      real(wp),intent(out)    :: WORK( * )
      end subroutine DTREXC
