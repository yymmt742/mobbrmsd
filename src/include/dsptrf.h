      pure subroutine DSPTRF( UPLO, N, AP, IPIV, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      integer,intent(in)      :: N
      integer,intent(out)     :: INFO
      integer,intent(out)     :: IPIV( * )
      real(wp),intent(inout)  :: AP( * )
      end subroutine DSPTRF
