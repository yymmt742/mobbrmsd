      pure subroutine DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: DIAG, UPLO
      integer,intent(in)      :: LDA, N
      integer,intent(out)     :: INFO
      real(wp),intent(inout)  :: A( LDA, * )
      end subroutine DTRTI2
