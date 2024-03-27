pure subroutine SSYMM(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL32
  character(*), intent(in) :: SIDE, UPLO
  real(rk), intent(in)     :: ALPHA, BETA
  integer, intent(in)      :: LDA, LDB, LDC, M, N
  real(rk), intent(in)     :: A(LDA, *), B(LDB, *)
  real(rk), intent(inout)  :: C(LDC, *)
end subroutine SSYMM

