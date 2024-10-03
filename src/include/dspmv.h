      pure subroutine DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: UPLO
      integer,intent(in)      :: INCX,INCY,N
      real(wp),intent(in)     :: ALPHA,BETA
      real(wp),intent(in)     :: AP(*),X(*)
      real(wp),intent(inout)  :: Y(*)
      end subroutine DSPMV
