      pure subroutine DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, &
     &                        INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      integer,intent(in)      :: LDA, LWORK, N
      integer,intent(out)     :: INFO
      real(wp),intent(inout)  :: A( LDA, * )
      real(wp),intent(out)    :: D( * ), E( * ), TAU( * ), WORK( * )
      end subroutine DSYTRD
