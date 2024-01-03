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
  integer(kind=ik), save     :: dmn = 0
  integer(kind=ik), save     :: nmem = 0
  integer(kind=ik), save     :: njob = 0
  logical, save              :: swap_y = .TRUE.
  type(symRMSD_input), save  :: inp
  type(symRMSD), allocatable :: sym(:)
!
contains
!
  subroutine add_molecule(m, n, s, sym)
    integer(kind=ik), intent(in) :: m
    integer(kind=ik), intent(in) :: n
    integer(kind=ik), intent(in) :: s
    integer(kind=ik), intent(in) :: sym(*)
    type(mol_block)         :: b_
    b_ = mol_block(1, s, m, n, m, n)
    call inp%add_molecule(b_, sym)
    dmn = inp%blk%mn * inp%blk%d
  end subroutine add_molecule
!
  subroutine setup()
    integer(kind=ik) :: i
    if (ALLOCATED(sym)) then
      deallocate (sym)
    end if
    !$omp parallel
    njob = MAX(omp_get_num_threads(), 1)
    !$omp end parallel
    ALLOCATE(sym(njob))
    do concurrent(i=1:njob)
      sym(i) = symRMSD(inp)
    end do
    nmem = sym(1)%nmem
  end subroutine setup
!
  subroutine run(x, y, n, res)
    real(kind=rk), intent(in)    :: x(*)
    real(kind=rk), intent(inout) :: y(*)
    integer(kind=ik), intent(in) :: n
    real(kind=rk), intent(out)   :: res(n)
    real(kind=rk)                :: w(nmem, njob)
    integer(kind=ik)             :: i
!
    !$omp parallel do
    do i = 1, n
      block
        integer(kind=ik) :: ijob, pnt
        pnt = (i - 1) * dmn + 1
        ijob = omp_get_thread_num() + 1
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
