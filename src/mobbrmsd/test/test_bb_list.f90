program main
  use mod_params, only: D, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_bb_block
  use mod_rotation
  use mod_bb_list
  use mod_unittest
  use mod_testutil
  implicit none
  type(unittest) :: u
!
  call u%init('test bb_list')
!
  call test0()
  call test1(4, 3, 1, [0])
  call test1(4, 2, 2, [3, 2, 1, 4])
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    type(bb_block)         :: blk(2)
    type(bb_list) :: b
    real(RK)               :: X(D, 8 * 3 + 4 * 5), Y(D, 8 * 3 + 4 * 5)
    integer(IK)            :: i, nmem
!
    blk(1) = bb_block(8, 3, sym=RESHAPE([2, 3, 4, 5, 6, 7, 8, 1], [8, 1]))
    blk(2) = bb_block(4, 5, sym=RESHAPE([1, 3, 2, 4, 4, 2, 3, 1], [4, 2]))
    b = bb_list(blk)
!
    X = sample(D, SIZE(X, 2))
    Y = X
!
    nmem = bb_list_memsize(b%q)
!
    block
      real(RK) :: w(nmem)
      do i = 1, 20
        call bb_list_setup(b%q, b%s, X, Y, w)
        call bb_list_run(b%q, b%s, w)
        print'(2F9.4,F9.1,2F9.4,8I3)', W(:4), EXP(W(4)), b%s(2:9)
        Y = 0.8 * Y + 0.2 * sample(D, SIZE(X, 2))
      end do
    end block
!
  end subroutine test0
!
  subroutine test1(m, n, s, sym)
    integer, intent(in)    :: m, n, s, sym(m * (s - 1))
    type(bb_block)         :: blk(2)
    type(bb_list) :: b
    real(RK), parameter    :: off = 5.0_RK
    real(RK)               :: X(D, m, n, 2), Y(D, m, n, 2)
    real(RK), allocatable  :: W(:)
    integer(IK)            :: i
!
    blk(1) = bb_block(m, n, sym=RESHAPE(sym, [m, s - 1]))
    blk(2) = bb_block(m, n, sym=RESHAPE(sym, [m, s - 1]))
    b = bb_list(blk)
!
    X = RESHAPE([sample(D, m * n) + off, sample(D, m * n) - off], SHAPE(X))
    Y = X
!
    allocate(W(bb_list_memsize(b%q)))
!
    do i = 1, 10
      call bb_list_setup(b%q, b%s, X, Y, W)
      call bb_list_run(b%q, b%s, w)
      call u%assert_almost_equal(W(1), brute_sd(m, n + n, s, sym, X, Y), 'minrmsd value')
      Y(:, :, :, 1) = 0.8 * Y(:, :, :, 1) + 0.2 * RESHAPE(sample(D, m * n) + off, [D, m, n])
      Y(:, :, :, 2) = 0.8 * Y(:, :, :, 2) + 0.2 * RESHAPE(sample(D, m * n) - off, [D, m, n])
    end do
!
  end subroutine test1
!
end program main

