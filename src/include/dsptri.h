      pure subroutine DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      integer,intent(in)      :: N
      integer,intent(out)     :: INFO
      integer,intent(in)      :: IPIV( * )
      real(wp),intent(inout)  :: AP( * )
      real(wp),intent(out)    :: WORK( * )
      end subroutine DSPTRI
