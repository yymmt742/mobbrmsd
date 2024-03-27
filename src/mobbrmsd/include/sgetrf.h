pure subroutine SGETRF(M, N, A, LDA, IPIV, INFO)
  use, intrinsic :: ISO_FORTRAN_ENV, only: rk => REAL32
  integer, intent(in)     :: LDA, M, N
  integer, intent(out)    :: INFO
  integer, intent(out)    :: IPIV(*)
  real(rk), intent(inout) :: A(LDA, *)
end subroutine SGETRF

