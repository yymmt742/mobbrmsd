program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_symRMSD
  use mod_unittest
  implicit none
  type(unittest) :: u
  integer, parameter :: NTEST=1
  integer            :: itest
!
  call u%init('test symRMSD')
  call test0()
!
  do itest = 1, NTEST
    call test1()
  end do
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    type(symRMSD_input)    :: inp
    integer(IK), parameter :: s(10) = [1, 3, 4, 5, 2, 2, 3, 5, 4, 1]
    call inp%add_molecule(mol_block(0, 1, 5, 4, 5, 4), s)
    call inp%add_molecule(mol_block(0, 2, 5, 4, 5, 4), s)
    call inp%add_molecule(mol_block(0, 3, 5, 4, 5, 4), s)
    print'(6i4)', inp%blk%b
    print*,size(inp%ms)
  end subroutine test0
!
  subroutine test1()
    integer, parameter         :: d = 3
    integer, parameter         :: s = 1
    integer, parameter         :: m = 5, n = 18, f = 12, g = 18
    type(mol_block), parameter :: b = mol_block(0, s, m, n, f, g)
    type(symRMSD_input)        :: inp
    type(symRMSD)              :: sr
    real(RK)                   :: X(d, m, n), Y(d, m, n)
    real(RK), allocatable      :: W(:)
    integer                    :: i, j, k
!
    inp%blk = mol_block_list(d, 1, [b])
!
    do k = -1, 1, 2
      do j = -1, 1
        do i = -1, 1
          X(:, :, (i + 1) * 6 + (j + 1) * 2 + (k + 1) / 2 + 1) = sample(d, m, [i, j, k])
          Y(:, :, (k + 1) / 2 * 9 + (j + 1) * 3 + (i + 1) + 1) = sample(d, m, [i, j, k])
        end do
      end do
    end do
!
    Y = 0.2 * X + 0.8 * Y
    sr = symRMSD(inp)
    allocate(w(sr%nmem))
    call sr%run(.true., X, Y, w)
    print *, SQRT(SUM((X - Y)**2) / (m * n)), sr%sd(W), sr%rmsd(W)
    print *, sr%search_ratio(W)
!
  end subroutine test1
!
  function sample(d, n, com) result(res)
    integer, intent(in)  :: d, n, com(d)
    real(RK)             :: res(d, n)
    integer              :: i
    call RANDOM_NUMBER(res)
    do concurrent(i=1:n)
      res(:, i) = res(:, i) + 6 * com
    enddo
  end function sample
!
end program main
