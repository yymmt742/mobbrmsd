      pure subroutine DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     :: K, LDA, LWORK, M, N
      integer,intent(out)    :: INFO
      real(wp),intent(in)    :: TAU( * )
      real(wp),intent(inout) :: A( LDA, * )
      real(wp),intent(out)   :: WORK( * )
      end subroutine DORGQL
