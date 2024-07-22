pure subroutine SORMQR(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL32
  character, intent(in)   :: SIDE, TRANS
  integer, intent(in)     :: M, N, K, LDA, LDC, LWORK
  integer, intent(out)    :: INFO
  real(rk), intent(in)    :: A(LDA, *), TAU(*)
  real(rk), intent(inout) :: C(LDC, *)
  real(rk), intent(out)   :: WORK(*)
end subroutine SORMQR

