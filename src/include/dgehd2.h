      pure subroutine DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     :: IHI, ILO, LDA, N
      integer,intent(out)    :: INFO
      real(wp),intent(inout) :: A( LDA, * )
      real(wp),intent(out)   :: TAU( * ), WORK( * )
      end subroutine DGEHD2
