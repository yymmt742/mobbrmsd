!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd
  use mod_dimspec_functions, only: D, setup_dimension
  use mod_params, only: &
 &      IK, &
 &      RK, &
 &      ONE => RONE, &
 &      ZERO => RZERO, &
 &      RHUGE
  use mod_bb_list
  use mod_bb_block
  use mod_mobbrmsd_header
  use mod_mobbrmsd_state
  implicit none
  public :: setup_dimension
  public :: mobbrmsd_input_init
  public :: mobbrmsd_input
  public :: mobbrmsd_input_add_molecule
  public :: mobbrmsd_input_destroy
  public :: mol_block_input
  public :: mol_block_input_add_molecule
  public :: mol_block_input_destroy
  public :: mobbrmsd
  public :: mobbrmsd_init
  public :: mobbrmsd_header
  public :: mobbrmsd_state
  public :: mobbrmsd_run
  public :: mobbrmsd_restart
  public :: mobbrmsd_is_finished
!
!| mol_block_input (for python interface)
  type mol_block_input
    sequence
    private
    integer(IK) :: m
    !! number of atoms per molecule
    integer(IK) :: n
    !! number of molecule
    integer(IK), allocatable :: sym(:, :)
    !! molecular symmetry
  end type mol_block_input
!
!| mobbrmsd_input
  type mobbrmsd_input
    sequence
    private
    type(mol_block_input), allocatable :: blk(:)
  end type mobbrmsd_input
!
!| mobbrmsd
  type mobbrmsd
    type(mobbrmsd_header) :: h
    !! mobbrmsd_header
    type(mobbrmsd_state)  :: s
    !! mobbrmsd_state
  contains
    final     :: mobbrmsd_destroy
  end type mobbrmsd
!
  interface mobbrmsd_init
    module procedure mobbrmsd_init_mobb, mobbrmsd_init_block
  end interface mobbrmsd_init
!
  interface mobbrmsd
    module procedure mobbrmsd_new, mobbrmsd_new_from_block
  end interface mobbrmsd
!
contains
!
!| init
  pure subroutine mobbrmsd_input_init(this)
    type(mobbrmsd_input), intent(inout) :: this
    !! self
    if (ALLOCATED(this%blk)) deallocate (this%blk)
  end subroutine mobbrmsd_input_init
!
!| add molecule
  pure subroutine mobbrmsd_input_add_molecule(this, m, n, sym)
    type(mobbrmsd_input), intent(inout) :: this
    !! self
    integer(IK), intent(in)             :: m
    !! number of atoms per molecule
    integer(IK), intent(in)             :: n
    !! number of molecule
    integer(IK), intent(in), optional   :: sym(:, :)
    !! molecular symmetry, sym(m, s-1)
    call mol_block_input_add_molecule(this%blk, m, n, sym)
  end subroutine mobbrmsd_input_add_molecule
!
!| destractor
  pure elemental subroutine mobbrmsd_input_destroy(this)
    type(mobbrmsd_input), intent(inout) :: this
    if (ALLOCATED(this%blk)) deallocate (this%blk)
  end subroutine mobbrmsd_input_destroy
!
! ------
!
!| add molecule
  pure subroutine mol_block_input_add_molecule(this, m, n, sym)
    type(mol_block_input), allocatable, intent(inout) :: this(:)
    !! mol_block_input array
    integer(IK), intent(in)            :: m
    !! number of atoms per molecule
    integer(IK), intent(in)            :: n
    !! number of molecule
    integer(IK), intent(in), optional  :: sym(:, :)
    !! molecular symmetry
    type(mol_block_input), allocatable :: blocks(:)
    integer(IK)                        :: nblock, i
!
    nblock = 1
    if (ALLOCATED(this)) nblock = nblock + SIZE(this)
!
    allocate (blocks(nblock))
    do concurrent(i=1:nblock - 1)
      blocks(i)%m = this(i)%m
      blocks(i)%n = this(i)%n
      call MOVE_ALLOC(from=this(i)%sym, to=blocks(i)%sym)
    end do
!
    blocks(nblock)%m = m
    blocks(nblock)%n = n
    if (PRESENT(sym)) then
      if (SIZE(sym) > 0) blocks(nblock)%sym = sym
    end if
!
    call MOVE_ALLOC(from=blocks, to=this)
!
  end subroutine mol_block_input_add_molecule
!
!| destractor
  pure elemental subroutine mol_block_input_destroy(this)
    type(mol_block_input), intent(inout) :: this
    if (ALLOCATED(this%sym)) deallocate (this%sym)
  end subroutine mol_block_input_destroy
!
! ------
!
!| constructor
  pure elemental function mobbrmsd_new(inp) result(res)
    type(mobbrmsd_input), intent(in) :: inp
    !! mobbrmsd_input
    type(mobbrmsd)                   :: res
    res = mobbrmsd_new_from_block(inp%blk)
  end function mobbrmsd_new
!
!| constructor, from mol_block_input. (for python interface)
  pure function mobbrmsd_new_from_block(blocks) result(res)
    type(mol_block_input), intent(in) :: blocks(:)
    !! mol_block_input array
    type(mobbrmsd)                    :: res
    call mobbrmsd_init(res, SIZE(blocks), blocks)
  end function mobbrmsd_new_from_block
!
!| constructor
  pure subroutine mobbrmsd_init_mobb(this, inp)
    type(mobbrmsd), intent(inout)    :: this
    !! self
    type(mobbrmsd_input), intent(in) :: inp
    !! mobbrmsd_input
    if (.not. ALLOCATED(inp%blk)) return
    call mobbrmsd_init_block(this, SIZE(inp%blk), inp%blk)
  end subroutine mobbrmsd_init_mobb
!
!| constructor
  pure subroutine mobbrmsd_init_block(this, nblock, blocks)
    type(mobbrmsd), intent(inout)     :: this
    !! self
    integer(IK), intent(in)           :: nblock
    !! mobbrmsd_input
    type(mol_block_input), intent(in) :: blocks(nblock)
    !! number of block
    type(bb_block)                    :: bbblk(nblock)
    !! blocks
    type(bb_list)                     :: bblst
    integer(IK)                       :: i
    do concurrent(i=1:nblock)
      if (ALLOCATED(blocks(i)%sym)) then
        call bb_block_init(bbblk(i), blocks(i)%m, blocks(i)%n, sym=blocks(i)%sym)
      else
        call bb_block_init(bbblk(i), blocks(i)%m, blocks(i)%n)
      end if
    end do
    call bb_list_init(bblst, nblock, bbblk)
    call mobbrmsd_header_init(this%h, SIZE(bblst%q), bblst%q, SIZE(bblst%s), bblst%s)
    call mobbrmsd_state_init(this%s, this%h)
    call bb_block_destroy(bbblk)
    call bb_list_destroy(bblst)
  end subroutine mobbrmsd_init_block
!
!| run mobbrmsd
  pure subroutine mobbrmsd_run( &
 &             header, &
 &             state, &
 &             X, Y, W, &
 &             cutoff, &
 &             difflim, &
 &             maxeval, &
 &             remove_com, &
 &             sort_by_g &
 &             )
    type(mobbrmsd_header), intent(in)   :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout) :: state
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(in)                :: X(*)
    !! reference coordinate
    real(RK), intent(in)                :: Y(*)
    !! target coordinate
    real(RK), intent(inout), optional   :: W(*)
    !! work array, must be > header%memsize()
    real(RK), intent(in), optional      :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional      :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional   :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional       :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional       :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
!
    if (PRESENT(W)) then
      call bb_list_setup(&
     &       header%q, &
     &       state%s,  &
     &       X,  &
     &       Y,  &
     &       W, &
     &       remove_com=remove_com, &
     &       sort_by_g=sort_by_g &
     &      )
      call mobbrmsd_restart( &
     &       header, &
     &       state,  &
     &       W, &
     &       cutoff=cutoff, &
     &       difflim=difflim, &
     &       maxeval=maxeval &
     &      )
    else
      block
        real(RK), allocatable :: T(:)
        allocate (T(header%memsize()))
        call bb_list_setup(&
       &       header%q, &
       &       state%s,  &
       &       X,  &
       &       Y,  &
       &       T, &
       &       remove_com=remove_com &
       &      )
        call mobbrmsd_restart( &
       &       header, &
       &       state,  &
       &       T, &
       &       cutoff=cutoff, &
       &       difflim=difflim, &
       &       maxeval=maxeval &
       &      )
      end block
    end if
  end subroutine mobbrmsd_run
!
!| run mobbrmsd
  pure subroutine mobbrmsd_restart( &
 &                  header, &
 &                  state, &
 &                  W, &
 &                  cutoff, &
 &                  difflim, &
 &                  maxeval &
 &                 )
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(inout)              :: W(*)
    !! work array, must be > header%memsize()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
    call bb_list_run( &
      &    header%q, &
      &    state%s, &
      &    W, &
      &    cutoff=cutoff, &
      &    difflim=difflim, &
      &    maxeval=maxeval &
      &   )
    call mobbrmsd_state_update(state, header, W)
  end subroutine mobbrmsd_restart
!
!| Returns bb process is finished.
  pure function mobbrmsd_is_finished(header, state) result(res)
    type(mobbrmsd_header), intent(in) :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(in)  :: state
    !! mobbrmsd_state
    logical                           :: res
    res = bb_list_is_finished(header%q, state%s)
  end function mobbrmsd_is_finished
!
  pure elemental subroutine mobbrmsd_destroy(this)
    type(mobbrmsd), intent(inout) :: this
  end subroutine mobbrmsd_destroy
end module mod_mobbrmsd

