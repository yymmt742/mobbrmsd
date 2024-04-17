!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd
!$ use omp_lib
  use blas_lapack_interface, only: D, setup_dimension
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, TEN => RTEN, LN_TO_L10, RHUGE
  use mod_iolib
  use mod_optarg
  use mod_bb_list
  use mod_bb_block
  implicit none
  public :: setup_dimension
  public :: mobbrmsd
  public :: mol_block_input
  public :: mol_block_input_add
  public :: mobbrmsd_n_dims
  public :: mobbrmsd_num_threads
  public :: mobbrmsd_header
  public :: mobbrmsd_state
  public :: mobbrmsd_run
  public :: mobbrmsd_restart
  public :: mobbrmsd_batch_run
  public :: mobbrmsd_nearest_neighbor
  public :: mobbrmsd_min_span_tree
  public :: mobbrmsd_swap_y
  public :: mobbrmsd_rotate_y
  public :: mobbrmsd_is_finished
!
  public :: INDEX_OF_RCP_NATOM
  public :: INDEX_OF_UPPERBOUND
  public :: INDEX_OF_LOWERBOUND
  public :: INDEX_OF_LOG_RATIO
  public :: INDEX_OF_N_EVAL
  public :: INDEX_TO_ROTMAT
!
  integer(IK), parameter :: INDEX_OF_RCP_NATOM  = 1
  !! Index of reciprocal natom of dumped state
  integer(IK), parameter :: INDEX_OF_UPPERBOUND = 2
  !! Index of upperbound of dumped state
  integer(IK), parameter :: INDEX_OF_LOWERBOUND = 3
  !! Index of lowerbound of dumped state
  integer(IK), parameter :: INDEX_OF_N_EVAL     = 4
  !! Index of n_eval of dumped state
  integer(IK), parameter :: INDEX_OF_LOG_RATIO  = 5
  !! Index of log_ratio of dumped state
  integer(IK), parameter :: INDEX_TO_ROTMAT     = 6
  !! Index to rotmatrix of dumped state
!
!| mol_block_input
  type mol_block_input
    private
    integer(IK) :: m
    !! number of atoms per molecule
    integer(IK) :: n
    !! number of molecule
    integer(IK), allocatable :: sym(:, :)
    !! molecular symmetry
  contains
    final :: mol_block_input_destroy
  end type mol_block_input
!
!| mobbrmsd_header
  type mobbrmsd_header
    private
    integer(IK)              :: d
    !! spatial dimension
    integer(IK), allocatable :: q(:)
    !! header array
    integer(IK), allocatable :: s(:)
    !! state template
  contains
    procedure :: n_dims       => mobbrmsd_header_n_dims
    procedure :: n_block      => mobbrmsd_header_n_block
    procedure :: n_atoms      => mobbrmsd_header_n_atoms
    procedure :: log_n_nodes  => mobbrmsd_header_log_n_nodes
    procedure :: frac_n_nodes => mobbrmsd_header_frac_n_nodes
    procedure :: exp_n_nodes  => mobbrmsd_header_exp_n_nodes
    procedure :: memsize      => mobbrmsd_header_memsize
    procedure :: dump         => mobbrmsd_header_dump
    procedure :: load         => mobbrmsd_header_load
    final     :: mobbrmsd_header_destroy
  end type mobbrmsd_header
!
!| mobbrmsd_state
  type mobbrmsd_state
    private
    integer(IK), allocatable :: s(:)
    real(RK), allocatable    :: z(:)
  contains
    procedure :: upperbound     => mobbrmsd_state_upperbound
    !! upperbound
    procedure :: lowerbound     => mobbrmsd_state_lowerbound
    !! lowerbound
    procedure :: sqrdev         => mobbrmsd_state_sqrdev
    !! rmsd
    procedure :: rmsd           => mobbrmsd_state_rmsd
    !! rmsd
    procedure :: n_eval         => mobbrmsd_state_n_eval
    !! number of lowerbound evaluation
    procedure :: eval_ratio     => mobbrmsd_state_eval_ratio
    !! ratio of evaluated node
    procedure :: log_eval_ratio => mobbrmsd_state_log_eval_ratio
    !! log ratio of evaluated node
    procedure :: dump           => mobbrmsd_state_dump
    !! dump current state
    procedure :: dump_real      => mobbrmsd_state_dump_real
    !! dump real part of current state
    procedure :: load           => mobbrmsd_state_load
    !! load state
    final     :: mobbrmsd_state_destroy
    !! destracter
  end type mobbrmsd_state
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
  interface mobbrmsd
    module procedure mobbrmsd_new_from_block
  end interface mobbrmsd
!
  interface
    include 'dgemm.h'
    include 'sgemm.h'
  end interface
!
contains
!
!| add molecule
  pure subroutine mol_block_input_add(this, m, n, sym)
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
      call move_alloc(from=this(i)%sym, to=blocks(i)%sym)
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
  end subroutine mol_block_input_add
!
!| destractor
  pure elemental subroutine mol_block_input_destroy(this)
    type(mol_block_input), intent(inout) :: this
    if (ALLOCATED(this%sym)) deallocate (this%sym)
  end subroutine mol_block_input_destroy
!
! ------
!
!| Returns spatial dimension
  pure elemental function mobbrmsd_header_n_dims(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    integer(IK)                        :: res
    res = this%d
  end function mobbrmsd_header_n_dims
!
!| Returns number of molecular blocks
  pure elemental function mobbrmsd_header_n_block(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    integer(IK)                        :: res
    res = bb_list_n_block(this%q)
  end function mobbrmsd_header_n_block
!
!| Returns header_memsize
  pure elemental function mobbrmsd_header_memsize(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    integer(IK)                        :: res
    res = bb_list_memsize(this%q)
  end function mobbrmsd_header_memsize
!
!| Returns n_atoms
  pure elemental function mobbrmsd_header_n_atoms(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    integer(IK)                        :: res
    res = bb_list_n_atoms(this%q)
  end function mobbrmsd_header_n_atoms
!
!| Returns log_n_nodes
  pure elemental function mobbrmsd_header_log_n_nodes(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    real(RK)                           :: res
    res = bb_list_log_n_nodes(this%q)
  end function mobbrmsd_header_log_n_nodes
!
!| returns number of nodes in fraction.
  pure function mobbrmsd_header_frac_n_nodes(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    real(RK)                           :: tmp, res
    tmp = LN_TO_L10 * bb_list_log_n_nodes(this%q)
    res = TEN**(tmp - REAL(INT(tmp), RK))
  end function mobbrmsd_header_frac_n_nodes
!
!| returns number of nodes in exp.
  pure function mobbrmsd_header_exp_n_nodes(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    integer(IK)                        :: res
    res = INT(LN_TO_L10 * bb_list_log_n_nodes(this%q), IK)
  end function mobbrmsd_header_exp_n_nodes
!
!| dump header as integer array
  pure function mobbrmsd_header_dump(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    integer(IK), allocatable           :: res(:)
    if (ALLOCATED(this%q) .and. ALLOCATED(this%s)) then
      res = [this%d, SIZE(this%q), this%q, this%s]
    else
      res = [this%d, 0]
    end if
  end function mobbrmsd_header_dump
!
!| load integer array as header
  pure subroutine mobbrmsd_header_load(this, q)
    class(mobbrmsd_header), intent(inout) :: this
    !! mobbrmsd_header
    integer(IK), intent(in)               :: q(:)
    !! header array
    integer(IK)                           :: sq, i
    this%d = q(1)
    sq = q(2)
    if (sq < 1) then
      this%q = [(0, i=1, 0)]
      this%s = [(0, i=1, 0)]
      return
    else
      this%q = q(3:2 + sq)
      this%s = q(3 + sq:)
    end if
  end subroutine mobbrmsd_header_load
!
!| destructer
  pure elemental subroutine mobbrmsd_header_destroy(this)
    type(mobbrmsd_header), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine mobbrmsd_header_destroy
!
! ------
!
!| returns upperbound
  pure elemental subroutine mobbrmsd_state_reset(this, header)
    class(mobbrmsd_state), intent(inout) :: this
    type(mobbrmsd_header), intent(in)    :: header
    real(RK)                             :: z(5 + header%n_dims()**2)
!
    this%s = header%s
!
    z(INDEX_OF_RCP_NATOM) = ONE / real(header%n_atoms(), RK)
    z(INDEX_OF_UPPERBOUND) = RHUGE
    z(INDEX_OF_LOWERBOUND) = ZERO
    z(INDEX_OF_LOG_RATIO) = -RHUGE
    z(INDEX_OF_N_EVAL) = ZERO
    call eye(header%n_dims(), z(INDEX_TO_ROTMAT))
!
    this%z = z
!
  contains
    pure subroutine eye(n_dims, e)
      integer(IK), intent(in) :: n_dims
      real(RK), intent(inout) :: e(n_dims, n_dims)
      integer(IK)             :: i, j
      do concurrent(i=1:n_dims, j=1:n_dims)
        e(i, j) = MERGE(ONE, ZERO, i == j)
      end do
    end subroutine eye
  end subroutine mobbrmsd_state_reset
!
!| returns upperbound
  pure elemental function mobbrmsd_state_upperbound(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = this%z(INDEX_OF_UPPERBOUND)
  end function mobbrmsd_state_upperbound
!
!| returns lowerbound
  pure elemental function mobbrmsd_state_lowerbound(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = this%z(INDEX_OF_LOWERBOUND)
  end function mobbrmsd_state_lowerbound
!
!| returns squared deviation
  pure elemental function mobbrmsd_state_sqrdev(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = MAX(ZERO, this%z(INDEX_OF_UPPERBOUND))
  end function mobbrmsd_state_sqrdev
!
!| returns rmsd
  pure elemental function mobbrmsd_state_rmsd(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = SQRT(this%z(INDEX_OF_RCP_NATOM) * MAX(ZERO, this%z(INDEX_OF_UPPERBOUND)))
  end function mobbrmsd_state_rmsd
!
!| returns n_eval
  pure elemental function mobbrmsd_state_n_eval(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    integer(IK)                       :: res
    res = NINT(this%z(INDEX_OF_N_EVAL), IK)
  end function mobbrmsd_state_n_eval
!
!| returns eval_ratio
  pure elemental function mobbrmsd_state_eval_ratio(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = EXP(this%z(INDEX_OF_LOG_RATIO))
  end function mobbrmsd_state_eval_ratio
!
!| returns log_eval_ratio
  pure elemental function mobbrmsd_state_log_eval_ratio(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = this%z(INDEX_OF_LOG_RATIO)
  end function mobbrmsd_state_log_eval_ratio
!
!| dump header as integer array
  pure function mobbrmsd_state_dump(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    integer(IK), allocatable          :: res(:)
    res = this%s
  end function mobbrmsd_state_dump
!
!| dump header as integer array
  pure function mobbrmsd_state_dump_real(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    real(RK), allocatable             :: res(:)
    res = this%z
  end function mobbrmsd_state_dump_real
!
!| load integer array as header
  pure subroutine mobbrmsd_state_load(this, ns, s, nz, z)
    class(mobbrmsd_state), intent(inout) :: this
    !! mobbrmsd_header
    integer(IK), intent(in)              :: ns, s(*), nz
    !! state integer array
    real(RK), intent(in)                 :: z(*)
    !! state real array
    this%s = s(:ns)
    this%z = z(:nz)
  end subroutine mobbrmsd_state_load
!
  pure elemental subroutine mobbrmsd_state_destroy(this)
    type(mobbrmsd_state), intent(inout) :: this
    if (ALLOCATED(this%s)) deallocate (this%s)
    if (ALLOCATED(this%z)) deallocate (this%z)
  end subroutine mobbrmsd_state_destroy
!
! ------
!
  pure function mobbrmsd_new_from_block(blocks) result(res)
    type(mol_block_input), intent(in) :: blocks(:)
    type(mobbrmsd)                    :: res
    integer(IK)                       :: nblock
!
    nblock = SIZE(blocks)
!
    block
      type(bb_block) :: bbblk(nblock)
      type(bb_list)  :: bblst
      integer(IK)    :: i
!
      do concurrent(i=1:nblock)
        if (ALLOCATED(blocks(i)%sym)) then
          bbblk(i) = bb_block(blocks(i)%m, blocks(i)%n, sym=blocks(i)%sym)
        else
          bbblk(i) = bb_block(blocks(i)%m, blocks(i)%n)
        end if
      end do
!
      bblst = bb_list(bbblk)
      res%h%d = D
      res%h%q = bblst%q
      res%h%s = bblst%s
      res%s%s = bblst%s
      call mobbrmsd_state_reset(res%s, res%h)
    end block
!
  end function mobbrmsd_new_from_block
!
!| returns spatial dimension
  pure elemental function mobbrmsd_n_dims() result(res)
    integer(IK) :: res
    res = D
  end function mobbrmsd_n_dims
!
  function mobbrmsd_num_threads() result(res)
    integer(IK)                        :: res
    !$omp parallel
    res = MAX(omp_get_num_threads(), 1)
    !$omp end parallel
  end function mobbrmsd_num_threads
!
!| run mobbrmsd
  pure subroutine mobbrmsd_run(header, state, &
 &                             X, Y, W, &
 &                             cutoff, difflim, maxeval)
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(in)                 :: X(*)
    !! reference coordinate
    real(RK), intent(in)                 :: Y(*)
    !! target coordinate
    real(RK), intent(inout), optional    :: W(*)
    !! work array, must be > header%memsize()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
!
    if (.not. ALLOCATED(header%q)) return
    if (.not. ALLOCATED(state%s)) return
    if (.not. ALLOCATED(state%z)) return
!
    if (PRESENT(W)) then
!
      call bb_list_setup(header%q, &
     &                   state%s,  &
     &                   X,  &
     &                   Y,  &
     &                   W)
      call mobbrmsd_restart(header, &
     &                      state,  &
     &                      W, &
     &                      cutoff=cutoff, &
     &                      difflim=difflim, &
     &                      maxeval=maxeval)
!
    else
!
      block
        real(RK), allocatable :: T(:)
        allocate (T(header%memsize()))
        call bb_list_setup(header%q, &
       &                   state%s,  &
       &                   X,  &
       &                   Y,  &
       &                   T)
        call mobbrmsd_restart(header, &
       &                      state,  &
       &                      T, &
       &                      cutoff=cutoff, &
       &                      difflim=difflim, &
       &                      maxeval=maxeval)
      end block
!
    end if
!
  end subroutine mobbrmsd_run
!
!| run mobbrmsd
  pure subroutine mobbrmsd_restart(header, state, W, &
 &                                 cutoff, difflim, maxeval)
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
!
    if (.not. ALLOCATED(header%q)) return
    if (.not. ALLOCATED(state%s)) return
    if (.not. ALLOCATED(state%z)) return
!
    call bb_list_run(header%q, state%s, W, cutoff=cutoff, difflim=difflim, maxeval=maxeval)
!
    state%z(INDEX_OF_UPPERBOUND) = W(bb_list_INDEX_OF_UPPERBOUND)
    state%z(INDEX_OF_LOWERBOUND) = W(bb_list_INDEX_OF_LOWERBOUND)
    state%z(INDEX_OF_N_EVAL)     = W(bb_list_INDEX_OF_N_EVAL)
    state%z(INDEX_OF_LOG_RATIO)  = LOG(W(bb_list_INDEX_OF_N_EVAL)) - W(bb_list_INDEX_TO_LOG_N_COMB)
!
    call bb_list_rotation_matrix(header%q, state%s, W, state%z(INDEX_TO_ROTMAT))
!
  end subroutine mobbrmsd_restart
!
  !| batch parallel run
  subroutine mobbrmsd_batch_run(n_target, header, state, &
  &                             X, Y, W, &
  &                             cutoff, difflim, maxeval, rotate_y)
    integer(IK), intent(in)              :: n_target
    !! number of target coordinates
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state(n_target)
    !! mobbrmsd_state, the result is contained in this structure.
    real(RK), intent(in)                 :: X(*)
    !! reference coordinate
    real(RK), intent(inout)              :: Y(*)
    !! target coordinate
    real(RK), intent(inout)              :: W(*)
   !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional        :: rotate_y
    !! The search ends when ncount exceeds maxiter.
    integer(kind=IK)                     :: i, itgt, ijob, ypnt, wpnt, ldy, ldw
!
    i = 0
    ldy = header%n_dims() * header%n_atoms()
    ldw = header%memsize()
!
    !$omp parallel private(itgt, ijob, ypnt, wpnt)
    do
      !$omp critical
      i = i + 1
      itgt = i
      !$omp end critical
      if (itgt > n_target) exit
      wpnt = ldw * omp_get_thread_num() + 1
      ypnt = (itgt - 1) * ldy + 1
      call mobbrmsd_run(header, state(itgt), X, Y(ypnt), W(wpnt), &
     &                  cutoff=cutoff, difflim=difflim, maxeval=maxeval)
      if (rotate_y) call mobbrmsd_rotate_y(header, state(itgt), Y(ypnt))
    end do
    !$omp end parallel
!
  end subroutine mobbrmsd_batch_run
!
  !| batch parallel nearest_neighbor
  subroutine mobbrmsd_nearest_neighbor(n_target, header, state, X, Y, W, &
 &                                     cutoff, difflim, maxeval, mask, nnval, nnidx)
    integer(IK), intent(in)             :: n_target
    !! number of target coordinates
    type(mobbrmsd_header), intent(in)   :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout) :: state(n_target)
    !! mobbrmsd_state, the result is contained in this structure.
    real(kind=RK), intent(in)           :: X(*)
   !! reference coordinate
    real(kind=RK), intent(in)           :: y(*)
   !! target coordinate
    real(kind=RK), intent(inout)        :: W(*)
   !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
    real(RK), intent(in), optional      :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional      :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional   :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional       :: mask(n_target)
    !! If .false., skip the calculation.
    real(RK), intent(out), optional     :: nnval
    !! nearest neighbor index
    integer(IK), intent(out), optional  :: nnidx
    !! nearest neighbor index
    real(kind=RK)                       :: cutoff_global, ub
    integer(kind=IK)                    :: i, itgt, ypnt, wpnt, ldy, ldw
!
    ldy = header%n_dims() * header%n_atoms()
    ldw = header%memsize()
!
    if(PRESENT(cutoff))then
      cutoff_global = MERGE(RHUGE, cutoff, cutoff < ZERO)
    else
      cutoff_global = RHUGE
    endif
!
    if (PRESENT(mask)) then
      cutoff_global = MIN(MINVAL(state%upperbound(), mask), cutoff_global)
    else
      cutoff_global = MIN(MINVAL(state%upperbound()), cutoff_global)
    end if
!
    i = 0
!
    !$omp parallel private(itgt, ub)
    do
      !$omp critical
      i = i + 1
      itgt = i
      ub = cutoff_global
      !$omp end critical
!
      if (itgt > n_target) exit
!
      if (PRESENT(mask)) then
        if (.not. mask(itgt)) cycle
      end if
      if (bb_list_is_finished(header%q, state(itgt)%s)) cycle
!
      if (ub < state(itgt)%lowerbound()) cycle
!
      wpnt = ldw * omp_get_thread_num() + 1
      ypnt = (itgt - 1) * ldy + 1
!
      call mobbrmsd_run(header, state(itgt), &
     &                  X, Y(ypnt), W(wpnt), &
     &                  cutoff=ub, &
     &                  difflim=difflim, &
     &                  maxeval=maxeval)
      ub = state(itgt)%upperbound()
!
      !$omp critical
      cutoff_global = MIN(ub, cutoff_global)
      !$omp end critical
    end do
    !$omp end parallel
!
    if (PRESENT(mask)) then
      if(PRESENT(nnval)) nnval = MINVAL(state%upperbound(), mask)
      if(PRESENT(nnidx)) nnidx = MINLOC(state%upperbound(), 1, mask)
    else
      if(PRESENT(nnval)) nnval = MINVAL(state%upperbound())
      if(PRESENT(nnidx)) nnidx = MINLOC(state%upperbound(), 1)
    end if
!
  end subroutine mobbrmsd_nearest_neighbor
!
  !| batch parallel nearest_neighbor
  subroutine mobbrmsd_min_span_tree(n_target, header, state, X, W, &
 &                                  cutoff, difflim, maxeval, &
 &                                  edges, weights, show_progress, unit &
 )
    integer(IK), intent(in)              :: n_target
    !! number of coordinates
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout)  :: state(n_target, n_target)
    !! mobbrmsd_state, the result is contained in this structure.
    real(kind=RK), intent(in)            :: X(*)
    !! coordinate sequence
    real(kind=RK), intent(inout)         :: W(*)
    !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
    real(RK), intent(in), optional       :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional       :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional    :: maxeval
    !! The search ends when ncount exceeds maxiter.
    integer(IK), intent(out), optional   :: edges(2, n_target - 1)
    !! minimum spanning tree edges
    real(RK), intent(out), optional      :: weights(n_target - 1)
    !! minimum spanning tree weights
    logical, intent(in), optional        :: show_progress
    !! if true, show progress bar
    integer(IK), intent(in), optional    :: unit
    !! device number for progress bar
    logical                              :: mask(n_target)
    integer(kind=IK)                     :: list(2, n_target)
    real(RK)                             :: vval(n_target - 1), nnval, cutoff_
    integer(kind=IK)                     :: i, j, xpnt, ldx, nnidx
!
    integer(kind=IK)                     :: unit_
    character(:), allocatable            :: deco, reset
!
    ldx = header%n_dims() * header%n_atoms()
!
    do concurrent(i=1:n_target, j=1:n_target)
      call mobbrmsd_state_reset(state(i, j), header)
    enddo
!
    mask(:) = .true.
    mask(1) = .false.
    list(1, 1) = 0
    list(2, 1) = 1
!
    if (PRESENT(unit)) then; unit_ = unit
    else; unit_ = mod_iolib_STDOUT
    end if
!
    deco = decorator(color='Y', carret=.true.)
    reset = mod_iolib_FS_RESET
!
    do j = 1, n_target - 1
!
      write (6, '(A, I8,A)', ADVANCE='NO') deco//'----------------', j, reset
      flush(6)
      vval(j) = RHUGE
      if (PRESENT(cutoff)) then
        cutoff_ = MIN(RHUGE, MAX(cutoff, ZERO))
      else
        cutoff_ = RHUGE
      end if
!
      do i = 1, j
!
        xpnt = (list(2, i) - 1) * ldx + 1
!
        call mobbrmsd_nearest_neighbor( &
       &  n_target, header, &
       &  state(:, list(2, i)), &
       &  X(xpnt), X, W, &
       &  cutoff=cutoff_, &
       &  difflim=difflim, &
       &  maxeval=maxeval, &
       &  mask=mask, &
       &  nnval=nnval,&
       &  nnidx=nnidx)
!
        if (nnval < vval(j)) then
          vval(j) = nnval
          cutoff_ = MIN(cutoff_, vval(j))
          list(1, j + 1) = list(2, i)
          list(2, j + 1) = nnidx
        end if
!
      end do
!
      mask(list(2, j + 1)) = .false.
!
    end do
!
    write (6, '(A)', ADVANCE='NO') decorate('', carret=.true., clear='l')
!
    do j = 1, n_target
      state(j, j)%z(INDEX_OF_UPPERBOUND) = ZERO
      do i = 1, j - 1
        if (state(i, j)%upperbound() < state(j, i)%upperbound())then
          state(j, i) = state(i, j)
        else
          state(i, j) = state(j, i)
        end if
      end do
    end do
!
    if(PRESENT(edges))then
      do concurrent(i=1:n_target - 1)
        edges(1, i) = list(1, i + 1)
        edges(2, i) = list(2, i + 1)
      enddo
    endif
!
    if(PRESENT(weights))then
      do concurrent(i=1:n_target - 1)
        weights(i) = vval(i)
      enddo
    endif
!
  end subroutine mobbrmsd_min_span_tree
!
!| swap target coordinate.
  pure subroutine mobbrmsd_swap_y(header, state, Y)
    class(mobbrmsd_header), intent(in) :: header
    !! mobbrmsd_header
    class(mobbrmsd_state), intent(in)  :: state
    !! mobbrmsd_state
    real(RK), intent(inout)            :: Y(*)
    !! target coordinate
!
    call bb_list_swap_y(header%q, state%s, Y(1))
!
  end subroutine mobbrmsd_swap_y
!
!| swap and rotate target coordinate.
  pure subroutine mobbrmsd_rotate_y(header, state, Y)
    class(mobbrmsd_header), intent(in) :: header
    !! mobbrmsd_header
    class(mobbrmsd_state), intent(in)  :: state
    !! mobbrmsd_state
    real(RK), intent(inout)            :: Y(*)
    !! target coordinate
!
    call bb_list_swap_y(header%q, state%s, Y)
    call rotate(header%n_dims(), header%n_atoms(), state%z(INDEX_TO_ROTMAT), Y)
!
  contains
!
    pure subroutine rotate(n_dims, n_atoms, R, Y)
      integer(IK), intent(in) :: n_dims, n_atoms
      real(RK), intent(in)    :: R(n_dims, n_dims)
      real(RK), intent(inout) :: Y(n_dims, n_atoms)
      real(RK)                :: T(n_dims, n_atoms)
      T = MATMUL(TRANSPOSE(R), Y)
      Y = T
    end subroutine rotate
!
  end subroutine mobbrmsd_rotate_y
!
!| Returns bb process is finished.
  pure function mobbrmsd_is_finished(header, state) result(res)
    class(mobbrmsd_header), intent(in) :: header
    !! mobbrmsd_header
    class(mobbrmsd_state), intent(in)  :: state
    !! mobbrmsd_state
    logical                            :: res
!
    res = bb_list_is_finished(header%q, state%s)
!
  end function mobbrmsd_is_finished
!
  pure elemental subroutine mobbrmsd_destroy(this)
    type(mobbrmsd), intent(inout) :: this
  end subroutine mobbrmsd_destroy
!
end module mod_mobbrmsd

