program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_symRMSD
  use mod_unittest
  implicit none
  type(unittest) :: u
! integer, parameter :: NTEST=25
! integer            :: itest
!
  call u%init('test symRMSD')
! do itest = 1, NTEST
    call test1()
! end do
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    integer, parameter         :: d = 3
    integer, parameter         :: s = 1
    integer, parameter         :: m = 50, n = 27, f = 12, g = 27
    type(mol_block), parameter :: b = mol_block(0, s, m, n, f, g)
    type(symRMSD_input)        :: inp
    type(symRMSD)              :: sr
    real(RK)                   :: X(d, m, n), Y(d, m, n), res
    integer                    :: i, j, k
!
    inp%blk = mol_block_list(d, 1, [b])
!
    do concurrent(k=-1:1, j=-1:1, i=-1:1)
      X(:, :, (i + 1) * 9 + (j + 1) * 3 + k + 2) = sample(d, m, [i, j, k])
      Y(:, :, (k + 1) * 9 + (j + 1) * 3 + i + 2) = sample(d, m, [i, j, k])
    end do
!
    sr = symRMSD(blk, inp)
    call sr%run(1, .true., X, Y, res)
!
  end subroutine test1
!
  pure function swp(d, m, n, per, sym, ms, X) result(res)
    integer(IK), intent(in)        :: d, m, n, per(:), sym(:)
    type(mol_symmetry), intent(in) :: ms
    real(RK), intent(in)           :: X(d, m, n)
    real(RK)                       :: tmp(d, m, n), res(d, m * n)
    integer(IK)                    :: i
    tmp = X
    do i = 1, SIZE(per)
      tmp(:, :, per(i)) = X(:, :, i)
      call ms%swap(d, tmp(:, :, per(i)), sym(i))
    end do
    res = RESHAPE(tmp, [d, m * n])
  end function swp
!
  pure function sd(d, X, Y) result(res)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: X(:, :), Y(:, :)
    real(RK)                :: C(d, d), R(d, d), W(100), res
    C = MATMUL(Y, TRANSPOSE(X))
    call estimate_rotation_matrix(d, SUM(X * X) + SUM(Y * Y), C, R, W)
    res = SUM(X**2) + SUM(Y**2) - 2 * SUM(C * R)
  end function sd
!
  function sample(d, n, com) result(res)
    integer, intent(in)  :: d, n, com(d)
    real(RK)             :: res(d, n)
    integer              :: i
    call RANDOM_NUMBER(res)
    do concurrent(i=1:n)
      res(:, i) = res(:, i) + 2 * com
    enddo
  end function sample
!
end program main
