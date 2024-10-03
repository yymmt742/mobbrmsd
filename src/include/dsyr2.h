      pure subroutine DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      real(wp),intent(in)     :: ALPHA
      integer,intent(in)      :: INCX,INCY,LDA,N
      real(wp),intent(in)     :: X(*),Y(*)
      real(wp),intent(inout)  :: A(LDA,*)
      end subroutine DSYR2
