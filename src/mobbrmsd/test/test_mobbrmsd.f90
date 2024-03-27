program main
  use mod_params, only: D, DD, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mobbrmsd
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer, parameter :: NTEST=1
  integer            :: itest
!
  call u%init('test symRMSD')
!
  call test0()
! call test2()
!
! do itest = 1, NTEST
!   call test1()
! end do
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    type(mobbrmsd_input)   :: inp
    integer(IK), parameter :: s(10) = [1, 3, 4, 5, 2, 2, 3, 5, 4, 1]
!   call inp%add_molecule(mol_block(0, 1, 5, 4, 4), s)
!   call inp%add_molecule(mol_block(0, 2, 5, 4, 4), s)
!   call inp%add_molecule(mol_block(0, 3, 5, 4, 4), s)
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

