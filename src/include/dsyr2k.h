      pure subroutine DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: TRANS,UPLO
      integer,intent(in)      :: K,LDA,LDB,LDC,N
      real(wp),intent(in)     :: ALPHA,BETA
      real(wp),intent(in)     :: A(LDA,*),B(LDB,*)
      real(wp),intent(inout)  :: C(LDC,*)
      end subroutine DSYR2K
