!| Module for manage D matrix.<br>
!  D(nx, ny) :: Residue matrices.<br>
!    - D_IJ = min_{R,s} Tr[C_IJs @ R]<br>
module mod_lowerbound
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_params, only: gemm, dot, copy
  use mod_mol_block
  use mod_rotation_matrix
  use mod_Hungarian
  implicit none
  private
  public :: worksize_lowerbound
  public :: lowerbound
!
contains
!
!| Inquire worksize of s_matrix.
  pure elemental function worksize_lowerbound(p, b) result(res)
    !| p :: level
    integer(IK), intent(in)    :: p
    !| b :: mol_block
    type(mol_block),intent(in) :: b
    integer(IK)                :: res
    integer(IK)                :: n1, n2
    n1 = MAX(b%x%n, b%y%n) - p
    n2 = MIN(b%x%n, b%y%n) - p
    res = MAX(worksize_Hungarian(n1, n2), 1 + worksize_sdmin())
  end function worksize_lowerbound
!
!| lowerbound function.
  pure subroutine lowerbound(p, b, G, C, D, W)
    !| p :: level
    integer(IK), intent(in)    :: p
    !| b :: mol_block
    type(mol_block),intent(in) :: b
    !| G :: partial variance, G
    real(RK), intent(in)       :: G
    !| C :: partial covariance matrix, C(d, d)
    real(RK), intent(in)       :: C(*)
    !| D :: residual matrix, D(n1, n2), here n1 = MAX(nx, ny) - p and n2 = MIN(nx, ny) - p.
    real(RK), intent(in)       :: D(*)
    !| W :: workarray
    real(RK), intent(inout)    :: W(*)
    integer(IK)                :: n1, n2
!
    W(1) = ZERO
    n1 = MAX(b%x%n, b%y%n) - p
    n2 = MIN(b%x%n, b%y%n) - p
!
    if (p < 0 .or. n2 < 0) return
    if (0 < n2) call Hungarian(n1, n2, D, W)
    if (0 < p) then
      call estimate_sdmin(G, C, W(2))
      W(1) = W(1) + W(2)
    end if
!
  end subroutine lowerbound
!
end module mod_lowerbound

