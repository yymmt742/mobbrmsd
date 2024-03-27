program main
  use mod_params, only: D, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mobbrmsd
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer, parameter :: NTEST=1
!
  call u%init('test symRMSD')
!
  call test0()
! call test2()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    type(mobbrmsd_input) :: inp
    type(mobbrmsd)       :: mobb
!
    inp = mobbrmsd_input()
    call inp%add_molecule(24, 1)
    call inp%add_molecule(4, 4, sym=RESHAPE([2, 1, 3, 4], [4, 1]))
    call inp%add_molecule(20, 3)
!
    mobb = mobbrmsd(inp)
!
    block
      type(mobbrmsd_header) :: h
      type(mobbrmsd_state)  :: s(3)
      real(RK)              :: X(D, 24 + 4 * 4 + 20 * 3)
      real(RK)              :: Y(D, 24 + 4 * 4 + 20 * 3)
      h = mobb%h
      s(:) = mobb%s
      X = sample(SIZE(X, 1), SIZE(X, 2))
      Y = 0.8 * X + 0.2 * sample(SIZE(Y, 1), SIZE(Y, 2))
      call mobbrmsd_run(h, s(1), X, Y)
    end block
!
  end subroutine test0
!
! subroutine test1()
!   integer, parameter         :: s = 1
!   integer, parameter         :: m = 5, n = 18, g = 18
!   type(mol_block), parameter :: b = mol_block(0, s, m, n, g)
!   type(symRMSD_input)        :: inp
!   type(symRMSD)              :: sr
!   real(RK)                   :: X(d, m, n), Y(d, m, n)
!   real(RK), allocatable      :: W(:)
!   integer                    :: i, j, k
!
!   inp%blk = mol_block_list(1, [b])
!
!   do k = -1, 1, 2
!     do j = -1, 1
!       do i = -1, 1
!         X(:, :, (i + 1) * 6 + (j + 1) * 2 + (k + 1) / 2 + 1) = sample(d, m, [i, j, k])
!         Y(:, :, (k + 1) / 2 * 9 + (j + 1) * 3 + (i + 1) + 1) = sample(d, m, [i, j, k])
!       end do
!     end do
!   end do
!
!   Y = 0.2 * X + 0.8 * Y
!   sr = symRMSD(inp)
!   allocate(w(sr%nmem))
!   call sr%run(.true., X, Y, w)
!   print *, SQRT(SUM((X - Y)**2) / (m * n)), sr%sd(W), sr%rmsd(W)
!   print *, sr%search_ratio(W)
!
! end subroutine test1
!
! subroutine test2()
!   integer, parameter         :: s = 1
!   integer, parameter         :: m = 5, n = 1, g = 18
!   type(mol_block), parameter :: b = mol_block(0, s, m, n, g)
!   type(symRMSD_input)        :: inp
!   type(symRMSD)              :: sr
!   real(RK)                   :: X(d, m, n), Y(d, m, n)
!   real(RK), allocatable      :: W(:)
!
!   inp%blk = mol_block_list(1, [b])
!
!   X(:, :, 1) = sample(d, m, [0, 0, 0])
!   Y(:, :, 1) = sample(d, m, [0, 0, 0])
!
!   Y = 0.2 * X + 0.8 * Y
!   sr = symRMSD(inp)
!   allocate(w(sr%nmem))
!   call sr%run(.true., X, Y, w)
!   print *, SQRT(SUM((X - Y)**2) / (m * n)), sr%sd(W), sr%rmsd(W)
!   print *, sr%search_ratio(W)
!
! end subroutine test2
!
end program main

