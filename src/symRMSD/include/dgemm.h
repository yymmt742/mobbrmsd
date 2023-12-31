pure subroutine DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  use mod_params, only :  wp=>R8
  real(wp),intent(in)     :: ALPHA,BETA
  integer,intent(in)      :: K,LDA,LDB,LDC,M,N
  character(*),intent(in) :: TRANSA,TRANSB
  real(wp),intent(in)     :: A(LDA,*),B(LDB,*)
  real(wp),intent(inout)  :: C(LDC,*)
end subroutine DGEMM
