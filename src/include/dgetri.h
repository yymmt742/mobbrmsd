      pure subroutine DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)      :: LDA, LWORK, N
      integer,intent(out)     :: INFO
      integer,intent(in)      :: IPIV( * )
      real(wp),intent(inout)  :: A( LDA, * )
      real(wp),intent(out)    :: WORK( * )
      end subroutine DGETRI
