pure subroutine DSYMM(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  use LA_CONSTANTS, only: wp => dp
  implicit none
  character(*), intent(in) :: SIDE, UPLO
  real(wp), intent(in)     :: ALPHA, BETA
  integer, intent(in)      :: LDA, LDB, LDC, M, N
  real(wp), intent(in)     :: A(LDA, *), B(LDB, *)
  real(wp), intent(inout)  :: C(LDC, *)
end subroutine DSYMM
