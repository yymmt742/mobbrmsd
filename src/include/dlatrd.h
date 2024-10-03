      pure subroutine DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      integer,intent(in)      :: LDA, LDW, N, NB
      real(wp),intent(inout)  :: A( LDA, * )
      real(wp),intent(out)    :: E( * ), TAU( * ), W( LDW, * )
      end subroutine DLATRD
