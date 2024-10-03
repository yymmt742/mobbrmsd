      pure subroutine DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      IMPLICIT NONE
      character(*),intent(in) :: JOBZ, UPLO
      integer,intent(in)      :: LDA, LWORK, N
      integer,intent(out)     :: INFO
      real(wp),intent(inout)  :: A( LDA, * )
      real(wp),intent(out)    :: W( * ), WORK( * )
      end subroutine DSYEV
