!| Module for manage D matrix.<br>
!  L(G, C, D) = SUM_{i=1,...,p} (G - 2tr[CR]) + min_{nu} SUM_{i=p+1,...,N} D_{i nu(p)}
module mod_lowerbound
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_params, only: gemm, dot, copy
  use mod_mol_block
  use mod_rotation_matrix
  use mod_Hungarian
  implicit none
  private
  public :: lowerbound_worksize
  public :: lowerbound
!
contains
!
!| Inquire worksize of lowerbound.
  pure elemental function lowerbound_worksize(b, p) result(res)
    type(mol_block),intent(in) :: b
    !! b :: mol_block
    integer(IK), intent(in)    :: p
    !! p :: level
    integer(IK)                :: res
    integer(IK)                :: n
    n = mol_block_nmol(b) - p
    res = MAX(Hungarian_worksize(n, n), 1 + sdmin_worksize())
  end function lowerbound_worksize
!
!| lowerbound function.
!  L(G, C, D) = SUM_{i=1,...,p} (G - 2tr[CR]) + min_{nu} SUM_{i=p+1,...,N} D_{i nu(p)}
  pure subroutine lowerbound(p, b, G, C, D, W)
    integer(IK), intent(in)    :: p
    !! p :: level
    type(mol_block),intent(in) :: b
    !! b :: mol_block
    real(RK), intent(in)       :: G
    !! G :: partial auto variance, G
    real(RK), intent(in)       :: C(*)
    !! C :: partial covariance matrix, C(d, d)
    real(RK), intent(in)       :: D(*)
    !! D :: residual matrix, D(n1, n2), here n1 = MAX(nx, ny) - p and n2 = MIN(nx, ny) - p.
    real(RK), intent(inout)    :: W(*)
    !! W :: workarray, must be SIZE(W) > lowerbound_worksize(p, b).
    integer(IK)                :: n
!
    W(1) = ZERO
    n = mol_block_nmol(b) - p
!
    if (p < 0 .or. n < 0) return
    if (0 < n) call Hungarian(n, n, D, W)
    if (0 < p) then
      call estimate_sdmin(G, C, W(2))
      W(1) = W(1) + W(2)
    end if
!
  end subroutine lowerbound
!
end module mod_lowerbound

