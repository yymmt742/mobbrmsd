module mod_symRMSD
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_branch_and_bound, only: branch_and_bound, DEF_maxeval, DEF_cutoff
  use mod_mol_block
  use mod_mol_symmetry
  implicit none
  public :: symRMSD
  public :: symRMSD_input
!
  type symRMSD_input
    integer(IK)                     :: maxeval = DEF_maxeval
    real(RK)                        :: cutoff  = DEF_cutoff
    type(mol_block_list)            :: blk
    type(mol_symmetry), allocatable :: ms(:)
  contains
    procedure :: add_molecule => symRMSD_input_add_molecule
    procedure :: clear        => symRMSD_input_clear
    final     :: symRMSD_input_destroy
  end type symRMSD_input
!
  type symRMSD
    private
    integer(IK), public    :: nmem = 0
    integer(IK), public    :: natm = 0
    type(branch_and_bound) :: bra
  contains
    procedure :: run             => symRMSD_run
    procedure :: sd              => symRMSD_sd
    procedure :: sd_with_error   => symRMSD_sd_with_error
    procedure :: rmsd            => symRMSD_rmsd
    procedure :: rmsd_with_error => symRMSD_rmsd_with_error
    procedure :: search_ratio    => symRMSD_search_ratio
    procedure :: clear           => symRMSD_clear
    final     :: symRMSD_destroy
  end type symRMSD
!
  interface symRMSD
    module procedure symRMSD_new
  end interface symRMSD
!
  interface
    include 'dcopy.h'
  end interface
!
contains
!
  pure subroutine symRMSD_input_add_molecule(this, b, s)
    class(symRMSD_input), intent(inout) :: this
    type(mol_block), intent(in)         :: b
    integer(IK), intent(in)             :: s(*)
    type(mol_symmetry), allocatable     :: ms(:)
    integer(IK)                         :: l, i, n, h(2)
!
    call this%blk%add_molecule(b)
!
    l = this%blk%nspecies()
    h(1) = this%blk%b(l)%m
    h(2) = this%blk%b(l)%s - 1
    n = h(1) * h(2)
!
    allocate (ms(l))
    do concurrent(i=1:l - 1)
      ms(i) = this%ms(i)
    end do
    ms(l) = mol_symmetry(RESHAPE(s(:n), h))
    call move_alloc(from=ms, to=this%ms)
!
  end subroutine symRMSD_input_add_molecule

  pure elemental subroutine symRMSD_input_clear(this)
    class(symRMSD_input), intent(inout) :: this
    call this%blk%clear()
    if (ALLOCATED(this%ms)) deallocate (this%ms)
  end subroutine symRMSD_input_clear
!
  pure elemental subroutine symRMSD_input_destroy(this)
    type(symRMSD_input), intent(inout) :: this
    call symRMSD_input_clear(this)
  end subroutine symRMSD_input_destroy
!
!!!
!
  pure function symRMSD_new(inp) result(res)
    type(symRMSD_input), intent(in) :: inp
    type(symRMSD)                   :: res
!
    if (ALLOCATED(inp%ms)) then
      res%bra = branch_and_bound(inp%blk, ms=inp%ms, maxeval=inp%maxeval, cutoff=inp%cutoff)
    else
      res%bra = branch_and_bound(inp%blk, maxeval=inp%maxeval, cutoff=inp%cutoff)
    end if
!
    res%nmem = res%bra%memsize
    res%natm = inp%blk%mn
!
  end function symRMSD_new
!
  pure subroutine symRMSD_run(this, swap_y, x, y, w)
    class(symRMSD), intent(in) :: this
    logical, intent(in)        :: swap_y
    real(RK), intent(in)       :: x(*)
    real(RK), intent(inout)    :: y(*)
    real(RK), intent(inout)    :: w(*)
!
    call this%bra%setup(x, y, w)
    call this%bra%run(w, swap_y)
    if (swap_y) call dcopy(this%bra%dmn, w(this%bra%yp), 1, y, 1)
!
  end subroutine symRMSD_run
!
  pure function symRMSD_sd(this, W) result(res)
    class(symRMSD), intent(in) :: this
    real(RK), intent(in)       :: w(*)
    real(RK)                   :: res
    res = W(this%bra%upperbound)
  end function symRMSD_sd
!
  pure function symRMSD_rmsd(this, W) result(res)
    class(symRMSD), intent(in) :: this
    real(RK), intent(in)       :: w(*)
    real(RK)                   :: res
    if (this%natm > 0) then
      res = SQRT(ABS(W(this%bra%upperbound)) / this%natm)
    else
      res = ZERO
    end if
  end function symRMSD_rmsd
!
  pure function symRMSD_sd_with_error(this, W) result(res)
    class(symRMSD), intent(in) :: this
    real(RK), intent(in)       :: w(*)
    real(RK)                   :: res(2)
    if (this%natm > 0) then
      res(1) = W(this%bra%lowerbound)
      res(2) = W(this%bra%upperbound)
    else
      res = ZERO
    end if
  end function symRMSD_sd_with_error
!
  pure function symRMSD_rmsd_with_error(this, W) result(res)
    class(symRMSD), intent(in) :: this
    real(RK), intent(in)       :: w(*)
    real(RK)                   :: res(2)
    if (this%natm > 0) then
      res(1) = SQRT(ABS(W(this%bra%lowerbound)) / this%natm)
      res(2) = SQRT(ABS(W(this%bra%upperbound)) / this%natm)
    else
      res = ZERO
    end if
  end function symRMSD_rmsd_with_error
!
  pure function symRMSD_search_ratio(this, W) result(res)
    class(symRMSD), intent(in) :: this
    real(RK), intent(in)       :: w(*)
    real(RK)                   :: res(3)
    res(1) = w(this%bra%ratio)
    res(2) = w(this%bra%lncmb)
    res(3) = w(this%bra%nsrch)
  end function symRMSD_search_ratio
!
  pure elemental subroutine symRMSD_clear(this)
    class(symRMSD), intent(inout) :: this
    this%nmem = 0
    call this%bra%clear()
  end subroutine symRMSD_clear
!
  pure elemental subroutine symRMSD_destroy(this)
    type(symRMSD), intent(inout) :: this
    call symRMSD_clear(this)
  end subroutine symRMSD_destroy
!
end module mod_symRMSD
