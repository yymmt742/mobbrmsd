pure subroutine DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
  use LA_CONSTANTS, only: wp => dp
  character(*), intent(in) :: TRANS
  integer, intent(in)      :: LDA, LDB, N, NRHS
  integer, intent(out)     :: INFO
  integer, intent(in)      :: IPIV(*)
  real(wp), intent(in)     :: A(LDA, *)
  real(wp), intent(inout)  :: B(LDB, *)
end subroutine DGETRS
