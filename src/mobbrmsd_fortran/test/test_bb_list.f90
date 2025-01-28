program main
  use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
  use mod_dimspec_functions, only: D, setup_dimension
  use mod_params, only: R8, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_bb_block
  use mod_rotation
  use mod_bb_list
  use mod_unittest
  use mod_testutil
  implicit none
  type(unittest) :: u
  character(32)  :: carg
  integer(IK)    :: d_
#ifdef USE_REAL32
  integer, parameter :: place = 1
#else
  integer, parameter :: place = 7
#endif
!
  call GET_COMMAND_ARGUMENT(1, value=carg)
  read (carg, *) d_
  call setup_dimension(d_)
!
  call u%init('test bb_list for (n,M,S)=(1,1,1)')
  call test1(1, 1, 1, [0])
  call u%init('test bb_list for (n,M,S)=(1,2,1)')
  call test1(1, 2, 1, [0])
  call u%init('test bb_list for (n,M,S)=(1,3,1)')
  call test1(1, 3, 1, [0])
  call u%init('test bb_list for (n,M,S)=(4,1,1)')
  call test1(4, 1, 1, [0])
  call u%init('test bb_list for (n,M,S)=(4,2,1)')
  call test1(4, 2, 1, [0])
  call u%init('test bb_list for (n,M,S)=(4,3,1)')
  call test1(4, 3, 1, [0])
  call u%init('test bb_list for (n,M,S)=(4,1,2)')
  call test1(4, 1, 2, [3, 2, 1, 4])
  call u%init('test bb_list for (n,M,S)=(4,2,2)')
  call test1(4, 2, 2, [3, 2, 1, 4])
  call u%init('test bb_list for (n,M,S)=(40,6,1)')
  call test1(40, 6, 1, [0])
!
  call u%init('test bb_list for {(n,M,S)}={(5,1,1), (5,4,1)}')
  call test2(5, 1, 1, [0], 5, 4, 1, [0])
  call u%init('test bb_list for {(n,M,S)}={(4,2,1), (5,2,1)}')
  call test2(4, 2, 1, [0], 5, 2, 1, [0])
  call u%init('test bb_list for {(n,M,S)}={(8,3,1), (4,1,2)}')
  call test2(8, 3, 1, [0], 4, 1, 1, [0])
  call u%init('test bb_list for {(n,M,S)}={(24,3,1), (24,4,1)}')
  call test2(24, 3, 1, [0], 24, 4, 1, [0])
!
  call u%init('test bb_list iterative for {(n,M,S)}={(5,1,1), (5,4,1)}')
  call test3(5, 1, 1, [0], 5, 4, 1, [0])
  call u%init('test bb_list iterative for {(n,M,S)}={(4,2,1), (5,2,1)}')
  call test3(4, 2, 1, [0], 5, 2, 1, [0])
  call u%init('test bb_list iterative for {(n,M,S)}={(8,3,1), (4,1,2)}')
  call test3(8, 3, 1, [0], 4, 1, 1, [0])
  call u%init('test bb_list iterative for {(n,M,S)}={(24,3,1), (24,4,1)}')
  call test3(24, 3, 1, [0], 24, 4, 1, [0])
!
  call u%finish_and_terminate(passing_score=0.95_R8)
!
contains
!
  subroutine test1(n, m, s, sym)
    integer, intent(in)   :: n, m, s, sym(n * (s - 1))
    type(bb_block)        :: blk(1)
    type(bb_list)         :: b
    real(RK)              :: X(D, n, m), Y(D, n, m), Z(D, n, m), R(D, D), sd, brute, rxz, nrm
    real(RK), allocatable :: W(:)
    integer(IK)           :: IX(n * m)
    integer(IK)           :: i
!
    call bb_block_init(blk(1), n, m, sym=RESHAPE(sym, [n, s - 1]))
    call bb_list_init(b, SIZE(blk), blk)
!
    X = sample(n, m)
    call centering(n, m, X)
    Y = X
    nrm = ONE / (n * m)
!
    allocate (W(bb_list_memsize(b%q)))
!
    do i = 1, 20
      Z = Y
      call bb_list_setup(b%q, b%s, X, Y, W)
      call bb_list_run(b%q, b%s, w)
      call bb_list_swap_y(b%q, b%s, Z)
      call bb_list_rotation_matrix(b%q, b%s, w, R)
      sd = w(bb_list_INDEX_TO_AUTOCORR) + w(bb_list_INDEX_TO_UPPERBOUND) * 2._RK
      brute = brute_sd(n, m, s, sym, X, Y)
      call centering(n, m, Z)
      rxz = SUM(([X] - [MATMUL(TRANSPOSE(R), RESHAPE(Z, [D, n * m]))])**2)
      call u%assert_almost_equal(sd * nrm, brute * nrm, 'minrmsd value', place=place)
      call u%assert_almost_equal(sd * nrm, rxz * nrm, 'sd vs rotmat', place=place)
      call bb_list_swap_indices(b%q, b%s, IX)
      call swap_YZ(n * m, ix, Y, Z)
      call centering(n, m, Z)
      rxz = SUM(([X] - [MATMUL(TRANSPOSE(R), RESHAPE(Z, [D, n * m]))])**2)
      call u%assert_almost_equal(sd * nrm, rxz * nrm, 'sd vs rotmat (by swap_indices)', place=place)
      Y = 0.5 * Y + 0.5 * sample(n, m)
    end do
!
  end subroutine test1
!
  subroutine test2(n1, m1, s1, sym1, n2, m2, s2, sym2)
    integer, intent(in)   :: n1, m1, s1, sym1(n1 * (s1 - 1))
    integer, intent(in)   :: n2, m2, s2, sym2(n2 * (s2 - 1))
    type(bb_block)        :: blk(2)
    type(bb_list)         :: b
    real(RK)              :: sd, brute, nrm
    real(RK)              :: X1(D, n1, m1), X2(D, n2, m2)
    real(RK)              :: Y1(D, n1, m1), Y2(D, n2, m2)
    real(RK)              :: X(D, n1 * m1 + n2 * m2)
    real(RK)              :: Z(D, n1 * m1 + n2 * m2)
    integer(IK)           :: IX(n1 * m1 + n2 * m2)
    integer(IK)           :: i
!
    call bb_block_init(blk(1), n1, m1, sym=RESHAPE(sym1, [n1, s1 - 1]))
    call bb_block_init(blk(2), n2, m2, sym=RESHAPE(sym2, [n2, s2 - 1]))
    call bb_list_init(b, SIZE(blk), blk)
!
    X1 = sample(n1, m1)
    X2 = sample(n2, m2)
    Y1 = X1
    Y2 = X2
    X = RESHAPE([X1, X2], SHAPE(X))
    call centering(SIZE(X, 2), X)
    nrm = ONE / (n1 * m1 + n2 * m2)
!
    block
      real(RK) :: W(bb_list_memsize(b%q)), R(D, D), rxz
      do i = 1, 50
        call bb_list_setup(b%q, b%s, [X1, X2], [Y1, Y2], W)
        call u%assert(.not. bb_list_is_finished(b%q, b%s), 'is not finished')
        call bb_list_run(b%q, b%s, W)
        call u%assert(bb_list_is_finished(b%q, b%s), 'is finished')
        sd = w(bb_list_INDEX_TO_AUTOCORR) + w(bb_list_INDEX_TO_UPPERBOUND) + w(bb_list_INDEX_TO_UPPERBOUND)
        brute = brute_sd_double(n1, m1, s1, sym1, n2, m2, s2, sym2, X1, Y1, X2, Y2)
        call u%assert_almost_equal(sd * nrm, brute * nrm, 'minrmsd value', place=place)
        Z = RESHAPE([Y1, Y2], SHAPE(Z))
        call centering(SIZE(Z, 2), Z)
        call bb_list_swap_y(b%q, b%s, Z)
        call bb_list_rotation_matrix(b%q, b%s, W, R)
        rxz = SUM((X - MATMUL(TRANSPOSE(R), Z))**2)
        call u%assert_almost_equal(sd * nrm, rxz * nrm, 'sd vs rotmat', place=place)
        call bb_list_swap_indices(b%q, b%s, IX)
        call swap_YZ(n1 * m1 + n2 * m2, ix, [Y1, Y2], Z)
        call centering(SIZE(Z, 2), Z)
        rxz = SUM((X - MATMUL(TRANSPOSE(R), Z))**2)
        call u%assert_almost_equal(sd * nrm, rxz * nrm, 'sd vs rotmat (by swap_indices)', place=place)
        Y1 = 0.5 * Y1 + 0.5 * sample(n1, m1)
        Y2 = 0.5 * Y2 + 0.5 * sample(n2, m2)
      end do
    end block
!
  end subroutine test2
!
  subroutine test3(n1, m1, s1, sym1, n2, m2, s2, sym2)
    integer, intent(in)   :: n1, m1, s1, sym1(n1 * (s1 - 1))
    integer, intent(in)   :: n2, m2, s2, sym2(n2 * (s2 - 1))
    type(bb_block)        :: blk(2)
    type(bb_list)         :: b
    real(RK)              :: sd, brute, nrm
    real(RK)              :: X1(D, n1, m1), X2(D, n2, m2)
    real(RK)              :: Y1(D, n1, m1), Y2(D, n2, m2)
    real(RK)              :: X(D, n1 * m1 + n2 * m2)
    real(RK)              :: Z(D, n1 * m1 + n2 * m2)
    integer(IK)           :: IX(n1 * m1 + n2 * m2)
    call bb_block_init(blk(1), n1, m1, sym=RESHAPE(sym1, [n1, s1 - 1]))
    call bb_block_init(blk(2), n2, m2, sym=RESHAPE(sym2, [n2, s2 - 1]))
    call bb_list_init(b, SIZE(blk), blk)
    X1 = sample(n1, m1)
    X2 = sample(n2, m2)
    Y1 = sample(n1, m1)
    Y2 = sample(n2, m2)
    X = RESHAPE([X1, X2], SHAPE(X))
    call centering(SIZE(X, 2), X)
    nrm = ONE / (n1 * m1 + n2 * m2)
    block
      real(RK) :: W(bb_list_memsize(b%q)), R(D, D), rxz
      call bb_list_setup(b%q, b%s, [X1, X2], [Y1, Y2], W)
      do while (.not. bb_list_is_finished(b%q, b%s))
        print '(A,3f9.4)', "# ", w(bb_list_INDEX_TO_N_EVAL), &
       &                         w(bb_list_INDEX_TO_UPPERBOUND), &
       &                         w(bb_list_INDEX_TO_LOWERBOUND)
        call bb_list_run(b%q, b%s, W, maxeval=0)
      end do
      sd = w(bb_list_INDEX_TO_AUTOCORR) + w(bb_list_INDEX_TO_UPPERBOUND) + w(bb_list_INDEX_TO_UPPERBOUND)
      brute = brute_sd_double(n1, m1, s1, sym1, n2, m2, s2, sym2, X1, Y1, X2, Y2)
      call u%assert_almost_equal(sd, brute, 'minrmsd value', place=place)
      Z = RESHAPE([Y1, Y2], SHAPE(Z))
      call centering(SIZE(Z, 2), Z)
      call bb_list_swap_y(b%q, b%s, Z)
      call bb_list_rotation_matrix(b%q, b%s, W, R)
      rxz = SUM((X - MATMUL(TRANSPOSE(R), Z))**2)
      call u%assert_almost_equal(sd * nrm, rxz * nrm, 'sd vs rotmat', place=place)
      call bb_list_swap_indices(b%q, b%s, IX)
      call swap_YZ(n1 * m1 + n2 * m2, ix, [Y1, Y2], Z)
      call centering(SIZE(Z, 2), Z)
      rxz = SUM((X - MATMUL(TRANSPOSE(R), Z))**2)
      call u%assert_almost_equal(sd * nrm, rxz * nrm, 'sd vs rotmat (by swap_indices)', place=place)
    end block
    FLUSH (OUTPUT_UNIT)
    FLUSH (ERROR_UNIT)
  end subroutine test3
!
  pure subroutine swap_YZ(n, ix, Y, Z)
    integer(IK), intent(in) :: n, ix(n)
    real(RK), intent(in)    :: Y(D, n)
    real(RK), intent(inout) :: Z(D, n)
    integer(IK)             :: i
    do concurrent(i=1:n)
      Z(:, i) = Y(:, ix(i))
    end do
  end subroutine swap_YZ
end program main

