!| tri_index
module mod_tril_index
  use ISO_FORTRAN_ENV, only: R4 => REAL32
  use mod_params, only: IK, RK
  implicit none
  private
  public :: cantor_pair, tril_index, trild_index
!
contains
!
  pure elemental subroutine cantor_pair(i, j, k)
    integer(IK), intent(in)    :: i, j
    integer(IK), intent(inout) :: k
    k = i + j
    k = k * (k - 3) / 2 + j + 1
  end subroutine cantor_pair
!
  pure elemental subroutine cantor_pair_inv(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    i = INT(0.5_R4 * (SQRT(real(8 * k - 7, R4)) - 1.0_R4), IK)
    j = k - i * (i + 1) / 2
    i = i - j + 2
  end subroutine cantor_pair_inv
!
  pure elemental subroutine tril_index(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    call cantor_pair_inv(k, i, j)
    j = i + j
    i = j - i
  end subroutine tril_index
!
  pure elemental subroutine trild_index(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    call cantor_pair_inv(k, i, j)
    j = i + j - 1
    i = j - i + 1
  end subroutine trild_index
!
end module mod_tril_index
