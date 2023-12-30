pure subroutine DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
  use mod_params, only :  wp=>R8
  character(*), intent(in) :: JOBU, JOBVT
  integer, intent(in)      :: LDA, LDU, LDVT, LWORK, M, N
  integer, intent(out)     :: INFO
  real(wp), intent(inout)  :: A(LDA, *)
  real(wp), intent(out)    :: S(*), U(LDU, *), VT(LDVT, *), WORK(*)
end subroutine DGESVD
