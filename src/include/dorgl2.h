      pure subroutine DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     :: K, LDA, M, N
      integer,intent(out)    :: INFO
      real(wp),intent(in)    :: TAU( * )
      real(wp),intent(inout) :: A( LDA, * )
      real(wp),intent(out)   :: WORK( * )
      end subroutine DORGL2
