!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd
  use blas_lapack_interface, only: D, setup_dimension
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, TEN => RTEN, LN_TO_L10, RHUGE
  use mod_bb_list
  use mod_bb_block
  implicit none
  public :: setup_dimension
  public :: mobbrmsd
  public :: mol_block_input
  public :: mol_block_input_add
  public :: mobbrmsd_n_dims
  public :: mobbrmsd_header
  public :: mobbrmsd_state
  public :: mobbrmsd_run
  public :: mobbrmsd_swap_y
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
    integer(IK), allocatable    :: q(:)
  contains
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
    real(RK)                 :: rcnatm
    real(RK)                 :: uppbou
    real(RK)                 :: lowbou
    real(RK)                 :: lograt
    real(RK)                 :: numevl
  contains
    procedure :: upperbound     => mobbrmsd_state_upperbound
    !! upperbound
    procedure :: lowerbound     => mobbrmsd_state_lowerbound
    !! lowerbound
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
contains
!
  pure subroutine mol_block_input_add(this, m, n, sym)
    type(mol_block_input), allocatable, intent(inout) :: this(:)
    !! this
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
    if (PRESENT(sym)) blocks(nblock)%sym = sym
!
    call MOVE_ALLOC(from=blocks, to=this)
!
  end subroutine mol_block_input_add
!
  pure elemental subroutine mol_block_input_destroy(this)
    type(mol_block_input), intent(inout) :: this
    if (ALLOCATED(this%sym)) deallocate (this%sym)
  end subroutine mol_block_input_destroy
!
! ------
!
!| Returns number of molecular blocks
  pure elemental function mobbrmsd_header_n_block(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    integer(IK)                        :: res
    res = bb_list_n_block(this%q)
  end function mobbrmsd_header_n_block
!
  pure elemental function mobbrmsd_header_memsize(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    integer(IK)                        :: res
    res = bb_list_memsize(this%q)
  end function mobbrmsd_header_memsize
!
  pure elemental function mobbrmsd_header_n_atoms(this) result(res)
    class(mobbrmsd_header), intent(in) :: this
    !! mobbrmsd_header
    integer(IK)                        :: res
    res = bb_list_n_atoms(this%q)
  end function mobbrmsd_header_n_atoms
!
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
    res = this%q
  end function mobbrmsd_header_dump
!
!| load integer array as header
  pure subroutine mobbrmsd_header_load(this, q)
    class(mobbrmsd_header), intent(inout) :: this
    !! mobbrmsd_header
    integer(IK), intent(in)               :: q(:)
    !! header array
    this%q = q
  end subroutine mobbrmsd_header_load
!
!| destructer
  pure elemental subroutine mobbrmsd_header_destroy(this)
    type(mobbrmsd_header), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine mobbrmsd_header_destroy
!
! ------
!
  pure elemental function mobbrmsd_state_upperbound(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = MAX(ZERO, this%uppbou)
  end function mobbrmsd_state_upperbound
!
  pure elemental function mobbrmsd_state_lowerbound(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = MAX(ZERO, this%lowbou)
  end function mobbrmsd_state_lowerbound
!
  pure elemental function mobbrmsd_state_rmsd(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = SQRT(this%rcnatm * MAX(ZERO, this%uppbou))
  end function mobbrmsd_state_rmsd
!
  pure elemental function mobbrmsd_state_n_eval(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    integer(IK)                       :: res
    res = NINT(this%numevl, IK)
  end function mobbrmsd_state_n_eval
!
  pure elemental function mobbrmsd_state_eval_ratio(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = EXP(this%lograt)
  end function mobbrmsd_state_eval_ratio
!
  pure elemental function mobbrmsd_state_log_eval_ratio(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = this%lograt
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
    real(RK)                          :: res(5)
    res = [this%rcnatm, this%uppbou, this%lowbou, this%lograt, this%numevl]
  end function mobbrmsd_state_dump_real
!
!| load integer array as header
  pure subroutine mobbrmsd_state_load(this, ns, s, z)
    class(mobbrmsd_state), intent(inout) :: this
    !! mobbrmsd_header
    integer(IK), intent(in)              :: ns, s(*)
    !! state integer array
    real(RK), intent(in)                 :: z(*)
    !! state real array
    this%s = s(:ns)
    this%rcnatm = z(1)
    this%uppbou = z(2)
    this%lowbou = z(3)
    this%lograt = z(4)
    this%numevl = z(5)
  end subroutine mobbrmsd_state_load
!
  pure elemental subroutine mobbrmsd_state_destroy(this)
    type(mobbrmsd_state), intent(inout) :: this
    if (ALLOCATED(this%s)) deallocate (this%s)
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
!
      res%h%q = bblst%q
!
      res%s%s = bblst%s
      res%s%rcnatm = ONE / real(bb_list_n_atoms(bblst%q), RK)
      res%s%uppbou = RHUGE
      res%s%lowbou = ZERO
      res%s%lograt = ZERO
      res%s%numevl = ZERO
    end block
!
  end function mobbrmsd_new_from_block
!
  pure elemental function mobbrmsd_n_dims() result(res)
    integer(IK) :: res
    res = D
  end function mobbrmsd_n_dims

!| run mobbrmsd
  pure subroutine mobbrmsd_run(header, state, X, Y, W, cutoff, difflim, maxeval, rot)
    class(mobbrmsd_header), intent(in)   :: header
    !! mobbrmsd_header
    class(mobbrmsd_state), intent(inout) :: state
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
    real(RK), intent(inout), optional    :: rot(*)
    !! rotation matrix, if needed.
!
    if (.not. ALLOCATED(header%q)) return
    if (.not. ALLOCATED(state%s)) return
!
    if (PRESENT(W)) then
!
      call bb_list_setup(header%q, state%s, X(1), Y(1), W(1))
      call bb_list_run(header%q, state%s, W(1), cutoff=cutoff, difflim=difflim, maxeval=maxeval)
      state%uppbou = W(1)
      state%lowbou = W(2)
      state%numevl = W(3)
      state%lograt = W(4)
!
      if(PRESENT(rot)) call bb_list_rotation_matrix(header%q, state%s, W(1), rot(1))
!
    else
!
      block
        real(RK), allocatable :: T(:)
        allocate (T(header%memsize()))
        call bb_list_setup(header%q, state%s, X(1), Y(1), T(1))
        call bb_list_run(header%q, state%s, T, cutoff=cutoff, difflim=difflim, maxeval=maxeval)
        state%uppbou = T(1)
        state%lowbou = T(2)
        state%numevl = T(3)
        state%lograt = T(4)
!
        if(PRESENT(rot)) call bb_list_rotation_matrix(header%q, state%s, T(1), rot(1))
!
      end block
!
    end if
!
  end subroutine mobbrmsd_run
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
  pure elemental subroutine mobbrmsd_destroy(this)
    type(mobbrmsd), intent(inout) :: this
  end subroutine mobbrmsd_destroy
!
end module mod_mobbrmsd

