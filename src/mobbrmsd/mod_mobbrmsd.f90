!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd
  use mod_dimspec_functions, only: D, setup_dimension
  use mod_params, only: &
 &      IK, &
 &      RK, &
 &      ONE => RONE, &
 &      ZERO => RZERO, &
 &      TEN => RTEN, &
 &      LN_TO_L10, &
 &      RHUGE
  use mod_bb_list
  use mod_bb_block
  use mod_mobbrmsd_state
  implicit none
  public :: setup_dimension
  public :: mobbrmsd_input_init
  public :: mobbrmsd_input
  public :: mobbrmsd_input_add_molecule
  public :: mobbrmsd_input_destroy
  public :: mobbrmsd_state
  public :: mobbrmsd
  public :: mobbrmsd_init
  public :: mobbrmsd_run
  public :: mobbrmsd_restart
  public :: mobbrmsd_is_finished
  public :: mobbrmsd_n_dims
  public :: mobbrmsd_n_block
  public :: mobbrmsd_n_atoms
  public :: mobbrmsd_log_n_nodes
  public :: mobbrmsd_frac_n_nodes
  public :: mobbrmsd_exp_n_nodes
  public :: mobbrmsd_attributes
  public :: mobbrmsd_memsize
  public :: mobbrmsd_swap_and_rotation
  public :: mobbrmsd_dump
  public :: mobbrmsd_load
  public :: mobbrmsd_destroy
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
    private
    sequence
    integer(IK)              :: d
    !! spatial dimension
    integer(IK), allocatable :: q(:)
    !! header array
    integer(IK), allocatable :: s(:)
    !! state template
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
    this%d = D
    this%q = bblst%q
    this%s = bblst%s
    call bb_block_destroy(bbblk)
    call bb_list_destroy(bblst)
  end subroutine mobbrmsd_init_block
!
!| run mobbrmsd
  pure subroutine mobbrmsd_run( &
 &                  this, &
 &                  state, &
 &                  X, Y, W, &
 &                  cutoff, &
 &                  ub_cutoff, &
 &                  difflim, &
 &                  maxeval, &
 &                  remove_com, &
 &                  sort_by_g, &
 &                  get_rotation, &
 &                  difflim_absolute &
 &                 )
    type(mobbrmsd), intent(in)          :: this
    !! mobbrmsd
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
    real(RK), intent(in), optional      :: ub_cutoff
    !! The search ends when upperbound is determined to be greater than to ub_cutoff.
    real(RK), intent(in), optional      :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional   :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional       :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional       :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
    logical, intent(in), optional       :: get_rotation
    !! if true, calculate rotation matrix. default [.false.]
    logical, intent(in), optional       :: difflim_absolute
    !! if true, use absolute difflim. default [.false.]
    integer(IK)                         :: d_
!
    if (PRESENT(get_rotation)) then
      d_ = MERGE(this%d, 0, get_rotation)
    else
      d_ = 0
    end if
    call mobbrmsd_state_init(state, d_, mobbrmsd_n_atoms(this), this%s) ! if d_<1, omit rotation matrix calculation.
!
    if (PRESENT(W)) then
      call bb_list_setup(&
     &       this%q, &
     &       state%s,  &
     &       X,  &
     &       Y,  &
     &       W, &
     &       remove_com=remove_com, &
     &       sort_by_g=sort_by_g &
     &      )
      call mobbrmsd_restart( &
     &       this, &
     &       state,  &
     &       W, &
     &       cutoff=cutoff, &
     &       ub_cutoff=ub_cutoff, &
     &       difflim=difflim, &
     &       maxeval=maxeval, &
     &       difflim_absolute=difflim_absolute &
     &      )
    else
      block
        real(RK), allocatable :: T(:)
        allocate (T(mobbrmsd_memsize(this)))
        call bb_list_setup(&
       &       this%q, &
       &       state%s,  &
       &       X,  &
       &       Y,  &
       &       T, &
       &       remove_com=remove_com &
       &      )
        call mobbrmsd_restart( &
       &       this, &
       &       state,  &
       &       T, &
       &       cutoff=cutoff, &
       &       ub_cutoff=ub_cutoff, &
       &       difflim=difflim, &
       &       maxeval=maxeval, &
       &       difflim_absolute=difflim_absolute &
       &      )
      end block
    end if
  end subroutine mobbrmsd_run
!
!| run mobbrmsd
  pure subroutine mobbrmsd_restart( &
 &                  this, &
 &                  state, &
 &                  W, &
 &                  cutoff, &
 &                  ub_cutoff, &
 &                  difflim, &
 &                  maxeval, &
 &                  difflim_absolute &
 &                 )
    type(mobbrmsd), intent(in)           :: this
    !! mobbrmsd
    type(mobbrmsd_state), intent(inout)  :: state
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(inout)              :: W(*)
    !! work array, must be > header%memsize()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: ub_cutoff
    !! The search ends when lowerbound is determined to be greater than to ub_cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional        :: difflim_absolute
    !! if true, use absolute difflim. default [.false.]
    call bb_list_run( &
      &    this%q, &
      &    state%s, &
      &    W, &
      &    cutoff=cutoff, &
      &    ub_cutoff=ub_cutoff, &
      &    difflim=difflim, &
      &    maxeval=maxeval, &
      &    difflim_absolute=difflim_absolute &
      &   )
    associate ( &
   &   ac => state%z(mobbrmsd_state_INDEX_TO_AUTOCORR), &
   &   ub => state%z(mobbrmsd_state_INDEX_TO_UPPERBOUND), &
   &   lb => state%z(mobbrmsd_state_INDEX_TO_LOWERBOUND), &
   &   ne => state%z(mobbrmsd_state_INDEX_TO_N_EVAL), &
   &   lr => state%z(mobbrmsd_state_INDEX_TO_LOG_RATIO), &
   &   rt => mobbrmsd_state_INDEX_TO_ROTMAT, &
   &   bbac => W(bb_list_INDEX_TO_AUTOCORR), &
   &   bbub => W(bb_list_INDEX_TO_UPPERBOUND), &
   &   bblb => W(bb_list_INDEX_TO_LOWERBOUND), &
   &   bbne => W(bb_list_INDEX_TO_N_EVAL), &
   &   bbln => W(bb_list_INDEX_TO_LOG_N_COMB) &
    )
      ac = bbac
      ub = bbub
      lb = bblb
      ne = bbne
      lr = LOG(bbne) - bbln
      if (mobbrmsd_state_has_rotation_matrix(state)) call bb_list_rotation_matrix(this%q, state%s, W, state%z(rt))
    end associate
  end subroutine mobbrmsd_restart
!
!| Returns bb process is finished.
  pure function mobbrmsd_is_finished(this, state) result(res)
    type(mobbrmsd), intent(in)        :: this
    !! mobbrmsd
    type(mobbrmsd_state), intent(in)  :: state
    !! mobbrmsd_state
    logical                           :: res
    res = bb_list_is_finished(this%q, state%s)
  end function mobbrmsd_is_finished
!
!| Returns spatial dimension
  pure elemental function mobbrmsd_n_dims(this) result(res)
    type(mobbrmsd), intent(in) :: this
    !! this
    integer(IK)                :: res
    res = this%d
  end function mobbrmsd_n_dims
!
!| Returns number of molecular blocks
  pure elemental function mobbrmsd_n_block(this) result(res)
    type(mobbrmsd), intent(in) :: this
    !! this
    integer(IK)                :: res
    res = bb_list_n_block(this%q)
  end function mobbrmsd_n_block
!
!| Returns header_memsize
  pure elemental function mobbrmsd_memsize(this) result(res)
    type(mobbrmsd), intent(in) :: this
    !! this
    integer(IK)                :: res
    res = bb_list_memsize(this%q)
  end function mobbrmsd_memsize
!
!| Returns n_atoms
  pure elemental function mobbrmsd_n_atoms(this) result(res)
    type(mobbrmsd), intent(in) :: this
    !! this
    integer(IK)                :: res
    res = bb_list_n_atoms(this%q)
  end function mobbrmsd_n_atoms
!
!| Returns log_n_nodes
  pure elemental function mobbrmsd_log_n_nodes(this) result(res)
    type(mobbrmsd), intent(in) :: this
    !! this
    real(RK)                   :: res
    res = bb_list_log_n_nodes(this%q)
  end function mobbrmsd_log_n_nodes
!
!| returns number of nodes in fraction.
  pure function mobbrmsd_frac_n_nodes(this) result(res)
    type(mobbrmsd), intent(in) :: this
    !! this
    real(RK)                   :: tmp, res
    tmp = LN_TO_L10 * bb_list_log_n_nodes(this%q)
    res = TEN**(tmp - real(INT(tmp), RK))
  end function mobbrmsd_frac_n_nodes
!
!| returns number of nodes in exp.
  pure function mobbrmsd_exp_n_nodes(this) result(res)
    type(mobbrmsd), intent(in) :: this
    !! this
    integer(IK)                :: res
    res = INT(LN_TO_L10 * bb_list_log_n_nodes(this%q), IK)
  end function mobbrmsd_exp_n_nodes
!
!| swap and rotation y
  pure subroutine mobbrmsd_swap_and_rotation(this, state, X)
    type(mobbrmsd), intent(in)       :: this
    !! mobbrmsd header
    type(mobbrmsd_state), intent(in) :: state
    !! this
    real(RK), intent(inout)           :: X(*)
    !! coordinate
    associate (rt => mobbrmsd_state_INDEX_TO_ROTMAT)
      call bb_list_swap_y(this%q, state%s, X)
      if (mobbrmsd_state_has_rotation_matrix(state)) call rotate(this%d, mobbrmsd_n_atoms(this), state%z(rt), X)
    end associate
  contains
    pure subroutine rotate(n_dims, n_atoms, R, X)
      integer(IK), intent(in) :: n_dims, n_atoms
      real(RK), intent(in)    :: R(n_dims, n_dims)
      real(RK), intent(inout) :: X(n_dims, n_atoms)
      real(RK)                :: T(n_dims, n_atoms)
      T = MATMUL(TRANSPOSE(R), X)
      X = T
    end subroutine rotate
  end subroutine mobbrmsd_swap_and_rotation
!
!| attributes
  pure subroutine mobbrmsd_attributes(this, n_dim, n_atom, n_mem, n_header, n_int, n_float)
    type(mobbrmsd), intent(in)         :: this
    !! this
    integer(IK), intent(out), optional :: n_dim
    !! n_dim
    integer(IK), intent(out), optional :: n_atom
    !! length of memory
    integer(IK), intent(out), optional :: n_mem
    !! length of memory
    integer(IK), intent(out), optional :: n_header
    !! length of header array q
    integer(IK), intent(out), optional :: n_int
    !! length of state array s
    integer(IK), intent(out), optional :: n_float
    !! length of state array z
    if (PRESENT(n_dim)) n_dim = this%d
    if (PRESENT(n_atom)) n_atom = mobbrmsd_n_atoms(this)
    if (PRESENT(n_mem)) n_mem = mobbrmsd_memsize(this)
    if (PRESENT(n_header)) n_header = 3 + SIZE(this%q) + SIZE(this%s)
    call mobbrmsd_state_attributes(this%d, this%s, n_int, n_float)
  end subroutine mobbrmsd_attributes
!
!| dump header as integer array
  pure function mobbrmsd_dump(this) result(res)
    type(mobbrmsd), intent(in) :: this
    !! this
    integer(IK), allocatable   :: res(:)
    allocate (res, source=[this%d, SIZE(this%q), this%q, SIZE(this%s), this%s])
  end function mobbrmsd_dump
!
!| load integer array as header
  pure subroutine mobbrmsd_load(this, q)
    type(mobbrmsd), intent(inout) :: this
    !! this
    integer(IK), intent(in)       :: q(:)
    !! header array
    integer(IK)                   :: sq, ss, i
    this%d = q(1)
    sq = q(2)
    if (sq < 1) then
      this%q = [(0, i=1, 0)]
    else
      this%q = q(3:2 + sq)
    end if
    ss = q(3 + sq)
    if (ss < 1) then
      this%s = [(0, i=1, 0)]
    else
      this%s = q(3 + sq:2 + sq + ss)
    end if
  end subroutine mobbrmsd_load
!
  pure elemental subroutine mobbrmsd_destroy(this)
    type(mobbrmsd), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine mobbrmsd_destroy
!
end module mod_mobbrmsd

