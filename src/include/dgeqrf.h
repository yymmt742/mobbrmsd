pure subroutine DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL64
  integer, intent(in)     :: M, N, LDA, LWORK
  integer, intent(out)    :: INFO
  real(rk), intent(inout) :: A(LDA, *)
  real(rk), intent(out)   :: TAU(*), WORK(*)
end subroutine DGEQRF

