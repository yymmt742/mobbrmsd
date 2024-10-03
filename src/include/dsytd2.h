      pure subroutine DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      integer,intent(in)      :: LDA, N
      integer,intent(out)     :: INFO
      real(wp),intent(inout)  :: A( LDA, * )
      real(wp),intent(out)    :: D( * ), E( * ), TAU( * )
      end subroutine DSYTD2
