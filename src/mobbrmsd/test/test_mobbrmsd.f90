program main
  use blas_lapack_interface, only : D
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mobbrmsd
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test mobbrmsd for (n,M,S)=(1,1,1)')
  call test1(1, 1, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(1,2,1)')
  call test1(1, 2, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(1,3,1)')
  call test1(1, 3, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(4,1,1)')
  call test1(4, 1, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(4,3,1)')
  call test1(4, 2, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(4,3,1)')
  call test1(4, 3, 1, [0])
  call u%init('test mobbrmsd for (n,M,S)=(4,1,2)')
  call test1(4, 1, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd for (n,M,S)=(4,2,2)')
  call test1(4, 2, 2, [3, 2, 1, 4])
  call u%init('test mobbrmsd for (n,M,S)=(40,6,1)')
  call test1(40, 6, 1, [0])
!
  call u%init('test bb_list for {(n,M,S)}={(5,1,1), (5,4,1)}')
  call test2(5, 1, 1, [0], 5, 4, 1, [0])
  call u%init('test bb_list for {(n,M,S)}={(4,2,1), (5,2,1)}')
  call test2(4, 2, 1, [0], 5, 2, 1, [0])
  call u%init('test bb_list for {(n,M,S)}={(8,2,1), (4,2,2)}')
  call test2(8, 2, 1, [0], 4, 2, 2, [3, 2, 1, 4])
  call u%init('test bb_list for {(n,M,S)}={(24,3,1), (24,4,1)}')
  call test2(24, 3, 1, [0], 24, 4, 1, [0])
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(n, m, s, sym)
    integer, intent(in)    :: n, m, s, sym(n * (s - 1))
    type(mobbrmsd)         :: mobb
    type(mol_block_input), allocatable :: inp(:)
    real(RK)               :: X(D, n, m), Y(D, n, m)
    real(RK), allocatable  :: W(:)
    integer(IK)            :: i
!
    call mol_block_input_add(inp, n, m, sym=RESHAPE(sym, [n, s - 1]))
    mobb = mobbrmsd(inp)
!
    X = sample(n, m)
    Y = X
!
    allocate(W(mobb%h%memsize()))
!
    do i = 1, 20
      call mobbrmsd_run(mobb%h, mobb%s, X, Y, W)
      call u%assert_almost_equal(w(1), brute_sd(n, m, s, sym, X, Y), 'minrmsd value')
      Y = 0.8 * Y + 0.2 * sample(n, m)
    end do
!
    deallocate(inp)
!
  end subroutine test1
!
  subroutine test2(n1, m1, s1, sym1, n2, m2, s2, sym2)
    integer, intent(in)   :: n1, m1, s1, sym1(n1 * (s1 - 1))
    integer, intent(in)   :: n2, m2, s2, sym2(n2 * (s2 - 1))
    type(mol_block_input), allocatable :: inp(:)
    type(mobbrmsd)        :: mobb
    real(RK)              :: brute
    real(RK)              :: X1(D, n1, m1), X2(D, n2, m2)
    real(RK)              :: Y1(D, n1, m1), Y2(D, n2, m2)
    real(RK)              :: Y(D, m1 * n1 + m2 * n2)
    real(RK), allocatable :: W(:)
    integer(IK)           :: i
!
    call mol_block_input_add(inp, n1, m1, sym=RESHAPE(sym1, [n1, s1 - 1]))
    call mol_block_input_add(inp, n2, m2, sym=RESHAPE(sym2, [n2, s2 - 1]))
    mobb = mobbrmsd(inp)
!
    X1 = sample(n1, m1)
    X2 = sample(n2, m2)
    Y1 = X1
    Y2 = X2
!
    allocate(W(mobb%h%memsize()))
!
    do i = 1, 10
      Y = RESHAPE([Y1, Y2], SHAPE(Y))
      call mobbrmsd_run(mobb%h, mobb%s, [X1, X2], Y, W)
      brute = brute_sd_double(n1, m1, s1, sym1, n2, m2, s2, sym2, X1, Y1, X2, Y2)
      call u%assert_almost_equal(w(1), brute, 'minrmsd value')
      Y1 = 0.8 * Y1 + 0.2 * sample(n1, m1)
      Y2 = 0.8 * Y2 + 0.2 * sample(n2, m2)
    end do
!
  end subroutine test2
!
end program main

