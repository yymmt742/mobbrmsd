!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd_state
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, TEN => RTEN, LN_TO_L10, RHUGE
  use mod_mobbrmsd_header
  use mod_bb_list
  implicit none
  public :: mobbrmsd_state
! public :: mobbrmsd_state_INDEX_OF_RCP_NATOM
! public :: mobbrmsd_state_INDEX_OF_UPPERBOUND
! public :: mobbrmsd_state_INDEX_OF_LOWERBOUND
! public :: mobbrmsd_state_INDEX_OF_LOG_RATIO
! public :: mobbrmsd_state_INDEX_OF_N_EVAL
! public :: mobbrmsd_state_INDEX_TO_ROTMAT
!&<
  integer(IK), parameter :: mobbrmsd_state_INDEX_OF_RCP_N_ATOMS = 1
  !! Index of reciprocal natom of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_OF_UPPERBOUND  = 2
  !! Index of upperbound of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_OF_LOWERBOUND  = 3
  !! Index of lowerbound of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_OF_N_EVAL      = 4
  !! Index of n_eval of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_OF_LOG_RATIO   = 5
  !! Index of log_ratio of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_ROTMAT      = 6
  !! Index to rotmatrix of dumped state
!&>
!| mobbrmsd_state
  type mobbrmsd_state
    private
    integer(IK), allocatable, public :: s(:)
    real(RK), allocatable    :: z(:)
  contains
    procedure :: update                 => mobbrmsd_state_update
    !! update state
    procedure :: upperbound             => mobbrmsd_state_upperbound
    !! upperbound
    procedure :: lowerbound             => mobbrmsd_state_lowerbound
    !! lowerbound
    procedure :: sqrared_deviation      => mobbrmsd_state_sqrared_deviation
    !! sqrared_deviation
    procedure :: mean_sqrared_deviation => mobbrmsd_state_mean_sqrared_deviation
    !! sqrared_deviation
    procedure :: rmsd                   => mobbrmsd_state_rmsd
    !! rmsd
    procedure :: n_eval                 => mobbrmsd_state_n_eval
    !! number of lowerbound evaluation
    procedure :: eval_ratio             => mobbrmsd_state_eval_ratio
    !! ratio of evaluated node
    procedure :: log_eval_ratio         => mobbrmsd_state_log_eval_ratio
    !! log ratio of evaluated node
    procedure :: dump                   => mobbrmsd_state_dump
    !! dump current state
    procedure :: dump_real              => mobbrmsd_state_dump_real
    !! dump real part of current state
    procedure :: load                   => mobbrmsd_state_load
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
    real(RK)                          :: z(5 + header%n_dims()**2)
    associate ( &
   &   RN => mobbrmsd_state_INDEX_OF_RCP_N_ATOMS, &
   &   UB => mobbrmsd_state_INDEX_OF_UPPERBOUND, &
   &   LB => mobbrmsd_state_INDEX_OF_LOWERBOUND, &
   &   NE => mobbrmsd_state_INDEX_OF_N_EVAL, &
   &   LR => mobbrmsd_state_INDEX_OF_LOG_RATIO, &
   &   RT => mobbrmsd_state_INDEX_TO_ROTMAT &
    )
      z(RN) = ONE / real(header%n_atoms(), RK)
      z(UB) = RHUGE
      z(LB) = ZERO
      z(NE) = -RHUGE
      z(LR) = ZERO
      call eye(header%n_dims(), z(RT))
    end associate
!
    allocate (res%s, source=header%state_template())
    allocate (res%z, source=z)
!
  contains
!
    pure subroutine eye(n_dims, e)
      integer(IK), intent(in) :: n_dims
      real(RK), intent(inout) :: e(n_dims, n_dims)
      integer(IK)             :: i, j
      do concurrent(i=1:n_dims, j=1:n_dims)
        e(i, j) = MERGE(ONE, ZERO, i == j)
      end do
    end subroutine eye
!
  end function mobbrmsd_state_new
!
!| update mobbrmsd_state
  pure subroutine mobbrmsd_state_update(this, header, W)
    class(mobbrmsd_state), intent(inout) :: this
    !! mobbrmsd header
    type(mobbrmsd_header), intent(in)    :: header
    !! mobbrmsd header
    real(RK), intent(in)                 :: W(*)
    !! mobbrmsd workarray
    associate( &
   &   UB => mobbrmsd_state_INDEX_OF_UPPERBOUND, &
   &   LB => mobbrmsd_state_INDEX_OF_LOWERBOUND, &
   &   NE => mobbrmsd_state_INDEX_OF_N_EVAL, &
   &   LR => mobbrmsd_state_INDEX_OF_LOG_RATIO, &
   &   RT => mobbrmsd_state_INDEX_TO_ROTMAT, &
   &   BBUB => bb_list_INDEX_OF_UPPERBOUND, &
   &   BBLB => bb_list_INDEX_OF_LOWERBOUND, &
   &   BBNE => bb_list_INDEX_OF_N_EVAL, &
   &   BBLN => bb_list_INDEX_TO_LOG_N_COMB &
    )
!
      this%z(UB) = W(BBUB)
      this%z(LB) = W(BBLB)
      this%z(NE) = W(BBNE)
      this%z(LR) = LOG(W(BBNE)) - W(BBLN)
!
      call bb_list_rotation_matrix(header%q, this%s, W, this%z(RT))
    end associate
  end subroutine mobbrmsd_state_update
!
!| returns upperbound
  pure elemental function mobbrmsd_state_upperbound(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = this%z(mobbrmsd_state_INDEX_OF_UPPERBOUND)
  end function mobbrmsd_state_upperbound
!
!| returns lowerbound
  pure elemental function mobbrmsd_state_lowerbound(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    associate (UB => mobbrmsd_state_INDEX_OF_UPPERBOUND)
      res = this%z(UB)
    end associate
  end function mobbrmsd_state_lowerbound
!
!| returns squared deviation
  pure elemental function mobbrmsd_state_sqrared_deviation(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    associate (UB => mobbrmsd_state_INDEX_OF_UPPERBOUND)
      res = MAX(ZERO, this%z(UB))
    end associate
  end function mobbrmsd_state_sqrared_deviation
!
!| returns squared deviation
  pure elemental function mobbrmsd_state_mean_sqrared_deviation(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    associate ( &
   &  RN => mobbrmsd_state_INDEX_OF_RCP_N_ATOMS, &
   &  UB => mobbrmsd_state_INDEX_OF_UPPERBOUND &
   &  )
      res = MAX(ZERO, this%z(UB)) * this%z(RN)
    end associate
  end function mobbrmsd_state_mean_sqrared_deviation
!
!| returns rmsd
  pure elemental function mobbrmsd_state_rmsd(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    associate ( &
   &  RN => mobbrmsd_state_INDEX_OF_RCP_N_ATOMS, &
   &  UB => mobbrmsd_state_INDEX_OF_UPPERBOUND &
   &  )
      res = SQRT(this%z(RN) * MAX(ZERO, this%z(UB)))
    end associate
  end function mobbrmsd_state_rmsd
!
!| returns n_eval
  pure elemental function mobbrmsd_state_n_eval(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    integer(IK)                       :: res
    associate (NE => mobbrmsd_state_INDEX_OF_N_EVAL)
      res = NINT(this%z(NE), IK)
    end associate
  end function mobbrmsd_state_n_eval
!
!| returns eval_ratio
  pure elemental function mobbrmsd_state_eval_ratio(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    associate (LR => mobbrmsd_state_INDEX_OF_LOG_RATIO)
      res = EXP(this%z(LR))
    end associate
  end function mobbrmsd_state_eval_ratio
!
!| returns log_eval_ratio
  pure elemental function mobbrmsd_state_log_eval_ratio(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    associate (LR => mobbrmsd_state_INDEX_OF_LOG_RATIO)
      res = this%z(LR)
    end associate
  end function mobbrmsd_state_log_eval_ratio
!
!| dump header as integer array (for python interface api)
  pure function mobbrmsd_state_dump(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    integer(IK), allocatable          :: res(:)
    ALLOCATE(res, source=this%s)
  end function mobbrmsd_state_dump
!
!| dump header as integer array (for python interface api)
  pure function mobbrmsd_state_dump_real(this) result(res)
    class(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    real(RK), allocatable             :: res(:)
    ALLOCATE(res, source=this%z)
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
    ALLOCATE(s_, source=s)
    ALLOCATE(z_, source=z)
    call move_alloc(from=s_, to=this%s)
    call move_alloc(from=z_, to=this%z)
  end subroutine mobbrmsd_state_load
!
  pure elemental subroutine mobbrmsd_state_destroy(this)
    type(mobbrmsd_state), intent(inout) :: this
    if (ALLOCATED(this%s)) deallocate (this%s)
    if (ALLOCATED(this%z)) deallocate (this%z)
  end subroutine mobbrmsd_state_destroy
!
end module mod_mobbrmsd_state

