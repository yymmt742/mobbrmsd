!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_params, only: copy
  use mod_bb_list
  use mod_bb_block
  implicit none
  public :: mobbrmsd
  public :: mobbrmsd_input
  public :: mobbrmsd_header
  public :: mobbrmsd_state
!
  integer(IK), parameter :: DEF_maxeval = -1
  real(RK), parameter    :: DEF_cutoff  = RHUGE
!
  type mol_block_input
    private
    integer(IK) :: m
    integer(IK) :: n
    integer(IK), allocatable :: sym(:, :)
  contains
    final :: mol_block_input_destroy
  end type mol_block_input
!
  type mobbrmsd_input
    private
    type(mol_block_input), allocatable :: blocks(:)
  contains
    procedure :: add_molecule => mobbrmsd_input_add_molecule
    final     :: mobbrmsd_input_destroy
  end type mobbrmsd_input
!
  type mobbrmsd_header
    private
    integer(IK), allocatable :: q(:)
  contains
    procedure :: n_block => mobbrmsd_header_n_block
    procedure :: memsize => mobbrmsd_header_memsize
    final :: mobbrmsd_header_destroy
  end type mobbrmsd_header
!
  type mobbrmsd_state
    private
    integer(IK), allocatable :: s(:)
  contains
    final :: mobbrmsd_state_destroy
  end type mobbrmsd_state
!
  type mobbrmsd
    type(mobbrmsd_header) :: h
    type(mobbrmsd_state)  :: s
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
  interface mobbrmsd_input
    module procedure mobbrmsd_input_new
  end interface mobbrmsd_input
!
  interface mobbrmsd
    module procedure mobbrmsd_new
  end interface mobbrmsd
!
contains
!
  pure elemental subroutine mol_block_input_destroy(this)
    type(mol_block_input), intent(inout) :: this
    if (ALLOCATED(this%sym)) deallocate (this%sym)
  end subroutine mol_block_input_destroy
!
! ------
!
  pure elemental function mobbrmsd_input_new() result(res)
    type(mobbrmsd_input) :: res
    allocate (res%blocks(0))
  end function mobbrmsd_input_new
!
  pure subroutine mobbrmsd_input_add_molecule(this, m, n, sym)
    class(mobbrmsd_input), intent(inout) :: this
    !! this
    integer(IK), intent(in)              :: m
    !! number of atoms per molecule
    integer(IK), intent(in)              :: n
    !! number of molecule
    integer(IK), intent(in), optional    :: sym(:, :)
    !! molecular symmetry
    type(mol_block_input), allocatable   :: blocks(:)
    integer(IK)                          :: nblock, i
!
    nblock = SIZE(this%blocks) + 1
!
    allocate (blocks(nblock))
    do concurrent(i=1:nblock - 1)
      blocks(i) = this%blocks(i)
    end do
!
    blocks(nblock)%m = m
    blocks(nblock)%n = n
    if (PRESENT(sym)) blocks(nblock)%sym = sym
!
    call MOVE_ALLOC(from=blocks, to=this%blocks)
!
  end subroutine mobbrmsd_input_add_molecule
!
  pure elemental subroutine mobbrmsd_input_destroy(this)
    type(mobbrmsd_input), intent(inout) :: this
    if (ALLOCATED(this%blocks)) deallocate (this%blocks)
  end subroutine mobbrmsd_input_destroy
!
! ------
!
!| Returns number of molecular blocks
  pure elemental function mobbrmsd_header_n_block(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    integer(IK)               :: res
    res = bb_list_n_block(this%q)
  end function mobbrmsd_header_n_block
!
  pure elemental function mobbrmsd_header_memsize(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    integer(IK)               :: res
    res = bb_list_memsize(this%q)
  end function mobbrmsd_header_memsize
!
  pure elemental subroutine mobbrmsd_header_destroy(this)
    type(mobbrmsd_header), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine mobbrmsd_header_destroy
!
! ------
!
  pure elemental subroutine mobbrmsd_state_destroy(this)
    type(mobbrmsd_state), intent(inout) :: this
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine mobbrmsd_state_destroy
!
! ------
!
  pure function mobbrmsd_new(inp) result(res)
    type(mobbrmsd_input), intent(in) :: inp
    type(mobbrmsd)                   :: res
    integer(IK)                      :: nblock
!
    nblock = SIZE(inp%blocks)
!
    block
      type(bb_block) :: bbblk(nblock)
      type(bb_list)  :: bblst
      integer(IK)    :: i
!
      do concurrent(i=1:nblock)
        if (ALLOCATED(inp%blocks(i)%sym)) then
          bbblk(i) = bb_block(inp%blocks(i)%m, inp%blocks(i)%n, sym=inp%blocks(i)%sym)
        else
          bbblk(i) = bb_block(inp%blocks(i)%m, inp%blocks(i)%n)
        end if
      end do
!
      bblst = bb_list(bbblk)
      res%h%q = bblst%q
      res%s%s = bblst%s
!
    end block
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
    !call this%bra%clear()
  end subroutine mobbrmsd_clear
!
  pure elemental subroutine mobbrmsd_destroy(this)
    type(mobbrmsd), intent(inout) :: this
    call mobbrmsd_clear(this)
  end subroutine mobbrmsd_destroy
!
end module mod_mobbrmsd

