      pure subroutine DSPR(UPLO,N,ALPHA,X,INCX,AP)
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      real(wp),intent(in)     :: ALPHA
      integer,intent(in)      :: INCX,N
      real(wp),intent(in)     :: X(*)
      real(wp),intent(inout)  :: AP(*)
      end subroutine DSPR
