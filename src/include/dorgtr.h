      pure subroutine DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      integer,intent(in)      :: LDA, LWORK, N
      integer,intent(out)     :: INFO
      real(wp),intent(in)     :: TAU( * )
      real(wp),intent(inout)  :: A( LDA, * )
      real(wp),intent(out)    :: WORK( * )
      end subroutine DORGTR
