      pure subroutine DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      INTEGER,intent(in)      :: INCX,INCY,LDA,N
      real(wp),intent(in)     :: ALPHA,BETA
      real(wp),intent(in)     :: A(LDA,*),X(*)
      real(wp),intent(inout)  :: Y(*)
      end subroutine DSYMV
