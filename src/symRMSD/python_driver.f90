module driver
  !$ use omp_lib
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_symRMSD
  use mod_mol_block
  use mod_mol_symmetry
  implicit none
  private
  public swap_y
  public add_molecule
  public setup
  public run
  public clear
!
  integer(IK), save          :: dmn = 0
  integer(IK), save          :: nmem = 0
  integer(IK), save          :: njob = 0
  logical, save              :: swap_y = .TRUE.
  type(symRMSD_input), save  :: inp
  type(symRMSD), allocatable :: sym(:)
!
contains
!
  subroutine add_molecule(m, n, s, sym)
    integer(IK), intent(in) :: m
    integer(IK), intent(in) :: n
    integer(IK), intent(in) :: s
    integer(IK), intent(in) :: sym(*)
    type(mol_block)         :: b_
    b_ = mol_block(1, s, m, n, m, n)
    call inp%add_molecule(b_, sym)
    dmn = inp%blk%mn * inp%blk%d
  end subroutine add_molecule
!
  subroutine setup()
    integer(IK) :: i
    call clear()
    !$omp parallel
    njob = omp_get_num_threads()
    !$omp end parallel
    print*,njob
    ALLOCATE(sym(njob))
    do concurrent(i=1:njob)
      sym(i) = symRMSD(inp)
    end do
    nmem = sym(1)%nmem
  end subroutine setup
!
  subroutine run(x, y, n, res)
    real(RK), intent(in)    :: x(*)
    real(RK), intent(inout) :: y(*)
    integer(IK), intent(in) :: n
    real(RK), intent(out)   :: res(n)
    real(RK)                :: w(nmem, njob)
    integer(IK)             :: i
!
    !$omp parallel do
    do i = 1, n
      block
        integer(IK) :: ijob, pnt
        pnt = (i - 1) * dmn + 1
        ijob = omp_get_thread_num() + 1
        print *, i, "Hello! N =", omp_get_num_threads(), " and I am ", omp_get_thread_num()
        call sym(ijob)%run(swap_y, x, y(pnt), w(1, ijob), res(i))
      end block
    end do
    !$omp end parallel do

  end subroutine run
!
  subroutine clear()
    call inp%clear()
    if (ALLOCATED(sym)) then
      call sym%clear()
      deallocate (sym)
    end if
    nmem = 0
    njob = 0
    dmn = 0
  end subroutine clear
!
end module driver
