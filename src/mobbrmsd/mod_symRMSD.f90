!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_params, only: copy
  use mod_branch_and_bound, only: branch_and_bound, DEF_maxeval, DEF_cutoff
  use mod_mol_block
  implicit none
  public :: mobbrmsd
  public :: mobbrmsd_input
!
  type mobbrmsd_input
    integer(IK)                     :: maxeval = DEF_maxeval
    real(RK)                        :: cutoff  = DEF_cutoff
!   type(mol_block_list)            :: blk
  contains
    procedure :: add_molecule => mobbrmsd_input_add_molecule
    procedure :: clear        => mobbrmsd_input_clear
    final     :: mobbrmsd_input_destroy
  end type mobbrmsd_input
!
  type mobbrmsd
    private
    integer(IK), public    :: nmem = 0
    integer(IK), public    :: natm = 0
!   type(branch_and_bound) :: bra
  contains
    procedure :: run             => mobbrmsd_run
!   procedure :: sd              => mobbrmsd_sd
!   procedure :: sd_with_error   => mobbrmsd_sd_with_error
!   procedure :: rmsd            => mobbrmsd_rmsd
!   procedure :: rmsd_with_error => mobbrmsd_rmsd_with_error
!   procedure :: search_ratio    => mobbrmsd_search_ratio
    procedure :: clear           => mobbrmsd_clear
    final     :: mobbrmsd_destroy
  end type mobbrmsd
!
  interface mobbrmsd
    module procedure mobbrmsd_new
  end interface mobbrmsd
!
contains
!
  pure subroutine mobbrmsd_input_add_molecule(this, b, s)
    class(mobbrmsd_input), intent(inout) :: this
    type(mol_block), intent(in)         :: b
    integer(IK), intent(in)             :: s(*)
    integer(IK)                         :: l, i, n, h(2)
!
!   call this%blk%add_molecule(b)
!
!   l = this%blk%nspecies()
!   h(1) = this%blk%b(l)%m
!   h(2) = this%blk%b(l)%s - 1
!   n = h(1) * h(2)
!
!   allocate (ms(l))
!   do concurrent(i=1:l - 1)
!     ms(i) = this%ms(i)
!   end do
!   ms(l) = mol_symmetry(RESHAPE(s(:n), h))
!   call move_alloc(from=ms, to=this%ms)
!
  end subroutine mobbrmsd_input_add_molecule

  pure elemental subroutine mobbrmsd_input_clear(this)
    class(mobbrmsd_input), intent(inout) :: this
  end subroutine mobbrmsd_input_clear
!
  pure elemental subroutine mobbrmsd_input_destroy(this)
    type(mobbrmsd_input), intent(inout) :: this
    call mobbrmsd_input_clear(this)
  end subroutine mobbrmsd_input_destroy
!
!!!
!
  pure function mobbrmsd_new(inp) result(res)
    type(mobbrmsd_input), intent(in) :: inp
    type(mobbrmsd)                   :: res
!
!   if (ALLOCATED(inp%ms)) then
!     res%bra = branch_and_bound(inp%blk, ms=inp%ms, maxeval=inp%maxeval, cutoff=inp%cutoff)
!   else
!     res%bra = branch_and_bound(inp%blk, maxeval=inp%maxeval, cutoff=inp%cutoff)
!   end if
!
!   res%nmem = res%bra%memsize
!   res%natm = inp%blk%mn
!
  end function mobbrmsd_new
!
  pure subroutine mobbrmsd_run(this, swap_y, x, y, w)
    class(mobbrmsd), intent(in) :: this
    logical, intent(in)        :: swap_y
    real(RK), intent(in)       :: x(*)
    real(RK), intent(inout)    :: y(*)
    real(RK), intent(inout)    :: w(*)
!
!   call this%bra%setup(x, y, w)
!   call this%bra%run(w, swap_y)
!   if (swap_y) call dcopy(this%bra%dmn, w(this%bra%yp), 1, y, 1)
!
  end subroutine mobbrmsd_run
!
! pure function mobbrmsd_sd(this, W) result(res)
!   class(mobbrmsd), intent(in) :: this
!   real(RK), intent(in)       :: w(*)
!   real(RK)                   :: res
!   res = W(this%bra%upperbound)
! end function mobbrmsd_sd
!
! pure function mobbrmsd_rmsd(this, W) result(res)
!   class(mobbrmsd), intent(in) :: this
!   real(RK), intent(in)       :: w(*)
!   real(RK)                   :: res
!   if (this%natm > 0) then
!     res = SQRT(ABS(W(this%bra%upperbound)) / this%natm)
!   else
!     res = ZERO
!   end if
! end function mobbrmsd_rmsd
!
! pure function mobbrmsd_sd_with_error(this, W) result(res)
!   class(mobbrmsd), intent(in) :: this
!   real(RK), intent(in)       :: w(*)
!   real(RK)                   :: res(2)
!   if (this%natm > 0) then
!     res(1) = W(this%bra%lowerbound)
!     res(2) = W(this%bra%upperbound)
!   else
!     res = ZERO
!   end if
! end function mobbrmsd_sd_with_error
!
! pure function mobbrmsd_rmsd_with_error(this, W) result(res)
!   class(mobbrmsd), intent(in) :: this
!   real(RK), intent(in)       :: w(*)
!   real(RK)                   :: res(2)
!   if (this%natm > 0) then
!     res(1) = SQRT(ABS(W(this%bra%lowerbound)) / this%natm)
!     res(2) = SQRT(ABS(W(this%bra%upperbound)) / this%natm)
!   else
!     res = ZERO
!   end if
! end function mobbrmsd_rmsd_with_error
!
! pure function mobbrmsd_search_ratio(this, W) result(res)
!   class(mobbrmsd), intent(in) :: this
!   real(RK), intent(in)       :: w(*)
!   real(RK)                   :: res(3)
!   res(1) = w(this%bra%ratio)
!   res(2) = w(this%bra%lncmb)
!   res(3) = w(this%bra%nsrch)
! end function mobbrmsd_search_ratio
!
  pure elemental subroutine mobbrmsd_clear(this)
    class(mobbrmsd), intent(inout) :: this
    this%nmem = 0
    !call this%bra%clear()
  end subroutine mobbrmsd_clear
!
  pure elemental subroutine mobbrmsd_destroy(this)
    type(mobbrmsd), intent(inout) :: this
    call mobbrmsd_clear(this)
  end subroutine mobbrmsd_destroy
!
end module mod_mobbrmsd

