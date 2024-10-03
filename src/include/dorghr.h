      pure subroutine DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     :: IHI, ILO, LDA, LWORK, N
      integer,intent(out)    :: INFO
      real(wp),intent(in)    :: TAU( * )
      real(wp),intent(inout) :: A( LDA, * )
      real(wp),intent(out)   :: WORK( * )
      end subroutine DORGHR
