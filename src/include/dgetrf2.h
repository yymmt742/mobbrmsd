      pure subroutine DGETRF2( M, N, A, LDA, IPIV, INFO )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)      :: LDA, M, N
      integer,intent(out)     :: INFO
      integer,intent(out)     :: IPIV( * )
      real(wp),intent(inout)  :: A( LDA, * )
      end subroutine DGETRF2
