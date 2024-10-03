      pure subroutine DLANSY( NORM, UPLO, N, A, LDA, RES, WORK )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: NORM, UPLO
      integer,intent(in)      :: LDA, N
      real(wp),intent(in)     :: A( LDA, * )
      real(wp),intent(out)    :: WORK( * ), RES
      end subroutine DLANSY
