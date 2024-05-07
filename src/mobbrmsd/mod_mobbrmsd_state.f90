!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd_state
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, TEN => RTEN, LN_TO_L10, RHUGE
  use mod_mobbrmsd_header
  use mod_bb_list
  implicit none
  public :: mobbrmsd_state
  public :: mobbrmsd_state_INDEX_TO_RCP_N_ATOMS
  public :: mobbrmsd_state_INDEX_TO_AUTOCORR
  public :: mobbrmsd_state_INDEX_TO_UPPERBOUND
  public :: mobbrmsd_state_INDEX_TO_LOWERBOUND
  public :: mobbrmsd_state_INDEX_TO_N_EVAL
  public :: mobbrmsd_state_INDEX_TO_LOG_RATIO
  public :: mobbrmsd_state_INDEX_TO_ROTMAT
!&<
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_RCP_N_ATOMS = 1
  !! Index of reciprocal natom of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_AUTOCORR    = 2
  !! Index to auto correlation
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_UPPERBOUND  = 3
  !! Index of upperbound of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_LOWERBOUND  = 4
  !! Index of lowerbound of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_N_EVAL      = 5
  !! Index of n_eval of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_LOG_RATIO   = 6
  !! Index of log_ratio of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_ROTMAT      = 7
  !! Index to rotmatrix of dumped state
!&>
!| mobbrmsd_state
  type mobbrmsd_state
    integer(IK), allocatable :: s(:)
    real(RK), allocatable    :: z(:)
  contains
    procedure :: upperbound => mobbrmsd_state_upperbound
    !! upperbound
    procedure :: lowerbound => mobbrmsd_state_lowerbound
    !! lowerbound
    procedure :: autovariance => mobbrmsd_state_autovariance
    !! lowerbound
    procedure :: squared_deviation => mobbrmsd_state_squared_deviation
    !! sqrared_deviation
    procedure :: mean_squared_deviation => mobbrmsd_state_mean_squared_deviation
    !! sqrared_deviation
    procedure :: rmsd => mobbrmsd_state_rmsd
    !! rmsd
    procedure :: lowerbound_as_rmsd => mobbrmsd_state_lowerbound_as_rmsd
    !! lowerbound_as_rmsd
    procedure :: n_eval => mobbrmsd_state_n_eval
    !! number of lowerbound evaluation
    procedure :: eval_ratio => mobbrmsd_state_eval_ratio
    !! ratio of evaluated node
    procedure :: log_eval_ratio => mobbrmsd_state_log_eval_ratio
    !! log ratio of evaluated node
    procedure :: rotation => mobbrmsd_state_rotation
    !! rotate given coordinate
    procedure :: dump => mobbrmsd_state_dump
    !! dump current state
    procedure :: dump_real => mobbrmsd_state_dump_real
    !! dump real part of current state
    procedure :: load => mobbrmsd_state_load
    !! load state
    final     :: mobbrmsd_state_destroy
    !! destracter
  end type mobbrmsd_state
!
  interface mobbrmsd_state
    module procedure mobbrmsd_state_new
  end interface mobbrmsd_state
!
contains
! ------
!
!| returns upperbound
  pure elemental function mobbrmsd_state_new(header) result(res)
    type(mobbrmsd_header), intent(in) :: header
    !! mobbrmsd header
    type(mobbrmsd_state)              :: res
    real(RK)                          :: z(6 + header%n_dims()**2)
    associate ( &
   &   RN => mobbrmsd_state_INDEX_TO_RCP_N_ATOMS, &
   &   AC => mobbrmsd_state_INDEX_TO_AUTOCORR, &
   &   UB => mobbrmsd_state_INDEX_TO_UPPERBOUND, &
   &   LB => mobbrmsd_state_INDEX_TO_LOWERBOUND, &
   &   NE => mobbrmsd_state_INDEX_TO_N_EVAL, &
   &   LR => mobbrmsd_state_INDEX_TO_LOG_RATIO, &
   &   RT => mobbrmsd_state_INDEX_TO_ROTMAT &
    )
      z(RN) = ONE / real(header%n_atoms(), RK)
      z(AC) = ZERO
      z(UB) = RHUGE
      z(LB) = -RHUGE
      z(NE) = -RHUGE
      z(LR) = ZERO
      call eye(header%n_dims(), z(RT))
!
      allocate (res%s, source=header%state_template())
      allocate (res%z, source=z)
    end associate
  contains
    pure subroutine eye(n_dims, e)
      integer(IK), intent(in) :: n_dims
      real(RK), intent(inout) :: e(n_dims, n_dims)
      integer(IK)             :: i, j
      do concurrent(i=1:n_dims, j=1:n_dims)
        e(i, j) = MERGE(ONE, ZERO, i == j)
      end do
    end subroutine eye
  end function mobbrmsd_state_new
!
!| returns upperbound
  pure elemental function mobbrmsd_state_upperbound(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    res = this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND)
  end function mobbrmsd_state_upperbound
!
!| returns lowerbound
  pure elemental function mobbrmsd_state_lowerbound(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    res = this%z(mobbrmsd_state_INDEX_TO_LOWERBOUND)
  end function mobbrmsd_state_lowerbound
!
!| returns autovariance
  pure elemental function mobbrmsd_state_autovariance(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    res = this%z(mobbrmsd_state_INDEX_TO_AUTOCORR)
  end function mobbrmsd_state_autovariance
!
!| returns squared deviation
  pure elemental function mobbrmsd_state_squared_deviation(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  ac => this%z(mobbrmsd_state_INDEX_TO_AUTOCORR), &
   &  ub => this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND) &
   &  )
      res = MAX(ZERO, ac + ub + ub)
    end associate
  end function mobbrmsd_state_squared_deviation
!
!| returns squared deviation
  pure elemental function mobbrmsd_state_mean_squared_deviation(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  rn => this%z(mobbrmsd_state_INDEX_TO_RCP_N_ATOMS), &
   &  ac => this%z(mobbrmsd_state_INDEX_TO_AUTOCORR), &
   &  ub => this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND) &
   &  )
      res = MAX(ZERO, (ac + ub + ub) * rn)
    end associate
  end function mobbrmsd_state_mean_squared_deviation
!
!| returns rmsd
  pure elemental function mobbrmsd_state_rmsd(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  rn => this%z(mobbrmsd_state_INDEX_TO_RCP_N_ATOMS), &
   &  ac => this%z(mobbrmsd_state_INDEX_TO_AUTOCORR), &
   &  ub => this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND) &
   &  )
      res = SQRT(MAX(ZERO, (ac + ub + ub) * rn))
    end associate
  end function mobbrmsd_state_rmsd
!
!| returns lowerbound as rmsd
  pure elemental function mobbrmsd_state_lowerbound_as_rmsd(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  rn => this%z(mobbrmsd_state_INDEX_TO_RCP_N_ATOMS), &
   &  ac => this%z(mobbrmsd_state_INDEX_TO_AUTOCORR), &
   &  lb => this%z(mobbrmsd_state_INDEX_TO_LOWERBOUND) &
   &  )
      res = SQRT(MAX(ZERO, (ac + lb + lb) * rn))
    end associate
  end function mobbrmsd_state_lowerbound_as_rmsd
!
!| returns n_eval
  pure elemental function mobbrmsd_state_n_eval(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    integer(IK)                       :: res
    associate (ne => this%z(mobbrmsd_state_INDEX_TO_N_EVAL))
      res = NINT(ne, IK)
    end associate
  end function mobbrmsd_state_n_eval
!
!| returns eval_ratio
  pure elemental function mobbrmsd_state_eval_ratio(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (LR => mobbrmsd_state_INDEX_TO_LOG_RATIO)
      res = EXP(this%z(LR))
    end associate
  end function mobbrmsd_state_eval_ratio
!
!| returns eval_ratio
  pure subroutine mobbrmsd_state_rotation(this, header, X)
    class(mobbrmsd_state), intent(in) :: this
    !! this
    type(mobbrmsd_header), intent(in) :: header
    !! mobbrmsd header
    real(RK), intent(inout)           :: X(*)
    !! coordinate
    associate (rt => mobbrmsd_state_INDEX_TO_ROTMAT)
      call bb_list_swap_y(header%q, this%s, X)
      call rotate(header%n_dims(), header%n_atoms(), this%z(rt), X)
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
  end subroutine mobbrmsd_state_rotation
!
!| returns log_eval_ratio
  pure elemental function mobbrmsd_state_log_eval_ratio(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = this%z(mobbrmsd_state_INDEX_TO_LOG_RATIO)
  end function mobbrmsd_state_log_eval_ratio
!
!| dump header as integer array (for python interface api)
  pure function mobbrmsd_state_dump(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    integer(IK), allocatable          :: res(:)
    allocate (res, source=this%s)
  end function mobbrmsd_state_dump
!
!| dump header as integer array (for python interface api)
  pure function mobbrmsd_state_dump_real(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    real(RK), allocatable             :: res(:)
    allocate (res, source=this%z)
  end function mobbrmsd_state_dump_real
!
!| load integer array as header (for python interface api)
  pure subroutine mobbrmsd_state_load(this, s, z)
    class(mobbrmsd_state), intent(inout) :: this
    !! mobbrmsd_header
    integer(IK), intent(in)              :: s(:)
    !! state integer array
    real(RK), intent(in)                 :: z(:)
    !! state real array
    integer(IK), allocatable             :: s_(:)
    real(RK), allocatable                :: z_(:)
    allocate (s_, source=s)
    allocate (z_, source=z)
    call MOVE_ALLOC(from=s_, to=this%s)
    call MOVE_ALLOC(from=z_, to=this%z)
  end subroutine mobbrmsd_state_load
!
  pure elemental subroutine mobbrmsd_state_destroy(this)
    type(mobbrmsd_state), intent(inout) :: this
    if (ALLOCATED(this%s)) deallocate (this%s)
    if (ALLOCATED(this%z)) deallocate (this%z)
  end subroutine mobbrmsd_state_destroy
end module mod_mobbrmsd_state

