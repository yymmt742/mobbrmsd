!| molecular orientation corrected RMSD with branch-and-bound.
module mod_mobbrmsd_state
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, TWO => RTWO, TEN => RTEN, LN_TO_L10, RHUGE
  use mod_bb_list
  implicit none
  private
  public :: mobbrmsd_state
  public :: mobbrmsd_state_init
  public :: mobbrmsd_state_copy
  public :: mobbrmsd_state_upperbound
  public :: mobbrmsd_state_lowerbound
  public :: mobbrmsd_state_autovariance
  public :: mobbrmsd_state_squared_deviation
  public :: mobbrmsd_state_mean_squared_deviation
  public :: mobbrmsd_state_rmsd
  public :: mobbrmsd_state_lowerbound_as_rmsd
  public :: mobbrmsd_state_n_eval
  public :: mobbrmsd_state_eval_ratio
  public :: mobbrmsd_state_log_eval_ratio
  public :: mobbrmsd_state_rotation_matrix
  public :: mobbrmsd_state_dump
  public :: mobbrmsd_state_dump_real
  public :: mobbrmsd_state_load
  public :: mobbrmsd_state_destroy
  public :: mobbrmsd_state_INDEX_TO_AUTOCORR
  public :: mobbrmsd_state_INDEX_TO_UPPERBOUND
  public :: mobbrmsd_state_INDEX_TO_LOWERBOUND
  public :: mobbrmsd_state_INDEX_TO_N_EVAL
  public :: mobbrmsd_state_INDEX_TO_LOG_RATIO
  public :: mobbrmsd_state_INDEX_TO_ROTMAT
!&<
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_AUTOCORR    = 1
  !! Index to auto correlation
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_UPPERBOUND  = 2
  !! Index of upperbound of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_LOWERBOUND  = 3
  !! Index of lowerbound of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_N_EVAL      = 4
  !! Index of n_eval of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_LOG_RATIO   = 5
  !! Index of log_ratio of dumped state
  integer(IK), parameter :: mobbrmsd_state_INDEX_TO_ROTMAT      = 6
  !! Index to rotmatrix of dumped state
!&>
!| mobbrmsd_state
  type mobbrmsd_state
    sequence
    integer(IK)              :: d, n
    integer(IK), allocatable :: s(:)
    real(RK), allocatable    :: z(:)
  end type mobbrmsd_state
!
contains
!
!| Init
  pure subroutine mobbrmsd_state_init(this, d, n, s)
    type(mobbrmsd_state), intent(inout) :: this
    !! Self
    integer(IK), intent(in)             :: d
    integer(IK), intent(in)             :: n
    integer(IK), intent(in)             :: s(:)
    real(RK), allocatable               :: z(:)
    associate ( &
   &   AC => mobbrmsd_state_INDEX_TO_AUTOCORR, &
   &   UB => mobbrmsd_state_INDEX_TO_UPPERBOUND, &
   &   LB => mobbrmsd_state_INDEX_TO_LOWERBOUND, &
   &   NE => mobbrmsd_state_INDEX_TO_N_EVAL, &
   &   LR => mobbrmsd_state_INDEX_TO_LOG_RATIO, &
   &   RT => mobbrmsd_state_INDEX_TO_ROTMAT &
    )
      this%d = d
      this%n = n
      this%s = s
      allocate (z(5 + this%d**2))
      z(AC) = ZERO
      z(UB) = RHUGE
      z(LB) = -RHUGE
      z(NE) = -RHUGE
      z(LR) = ZERO
      call eye(this%d, z(RT))
      call MOVE_ALLOC(from=z, to=this%z)
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
  end subroutine mobbrmsd_state_init
!
!| Pure elemental copy, for NVHPC.
  pure elemental subroutine mobbrmsd_state_copy(lhs, rhs)
    type(mobbrmsd_state), intent(inout) :: lhs
    !! Self
    type(mobbrmsd_state), intent(in)    :: rhs
    if (ALLOCATED(rhs%s)) lhs%s = rhs%s
    if (ALLOCATED(rhs%z)) lhs%z = rhs%z
  end subroutine mobbrmsd_state_copy
!
!| returns upperbound
  pure elemental function mobbrmsd_state_upperbound(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    res = this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND)
  end function mobbrmsd_state_upperbound
!
!| returns lowerbound
  pure elemental function mobbrmsd_state_lowerbound(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    res = this%z(mobbrmsd_state_INDEX_TO_LOWERBOUND)
  end function mobbrmsd_state_lowerbound
!
!| returns autovariance
  pure elemental function mobbrmsd_state_autovariance(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    res = this%z(mobbrmsd_state_INDEX_TO_AUTOCORR)
  end function mobbrmsd_state_autovariance
!
!| returns squared deviation
  pure elemental function mobbrmsd_state_squared_deviation(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  ac => this%z(mobbrmsd_state_INDEX_TO_AUTOCORR), &
   &  ub => this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND) &
   &  )
      res = MAX(ZERO, ac + TWO * ub)
    end associate
  end function mobbrmsd_state_squared_deviation
!
!| returns squared deviation
  pure elemental function mobbrmsd_state_mean_squared_deviation(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  ac => this%z(mobbrmsd_state_INDEX_TO_AUTOCORR), &
   &  ub => this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND) &
   &  )
      res = ONE / real(bb_list_n_atoms(this%s), RK)
      res = MAX(ZERO, (ac + TWO * ub) * res)
    end associate
  end function mobbrmsd_state_mean_squared_deviation
!
!| returns rmsd
  pure elemental function mobbrmsd_state_rmsd(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  ac => this%z(mobbrmsd_state_INDEX_TO_AUTOCORR), &
   &  ub => this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND) &
   &  )
      res = ONE / real(bb_list_n_atoms(this%s), RK)
      res = SQRT(MAX(ZERO, (ac + TWO * ub) * res))
    end associate
  end function mobbrmsd_state_rmsd
!
!| returns lowerbound as rmsd
  pure elemental function mobbrmsd_state_lowerbound_as_rmsd(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  ac => this%z(mobbrmsd_state_INDEX_TO_AUTOCORR), &
   &  lb => this%z(mobbrmsd_state_INDEX_TO_LOWERBOUND) &
   &  )
      res = ONE / real(bb_list_n_atoms(this%s), RK)
      res = SQRT(MAX(ZERO, (ac + TWO * lb) * res))
    end associate
  end function mobbrmsd_state_lowerbound_as_rmsd
!
!| returns n_eval
  pure elemental function mobbrmsd_state_n_eval(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    integer(IK)                       :: res
    associate (ne => this%z(mobbrmsd_state_INDEX_TO_N_EVAL))
      res = NINT(ne, IK)
    end associate
  end function mobbrmsd_state_n_eval
!
!| returns eval_ratio
  pure elemental function mobbrmsd_state_eval_ratio(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (LR => mobbrmsd_state_INDEX_TO_LOG_RATIO)
      res = EXP(this%z(LR))
    end associate
  end function mobbrmsd_state_eval_ratio
!
!| rotation matrix
  pure subroutine mobbrmsd_state_rotation_matrix(this, R)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK), intent(inout)           :: R(*)
    !! coordinate
    integer                           :: n
    associate (rt => mobbrmsd_state_INDEX_TO_ROTMAT)
      n = SIZE(this%z(rt:))
      R(:n) = this%z(rt:)
    end associate
  end subroutine mobbrmsd_state_rotation_matrix
!
!| returns log_eval_ratio
  pure elemental function mobbrmsd_state_log_eval_ratio(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    real(RK)                          :: res
    res = this%z(mobbrmsd_state_INDEX_TO_LOG_RATIO)
  end function mobbrmsd_state_log_eval_ratio
!
!| dump header as integer array (for python interface api)
  pure function mobbrmsd_state_dump(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    integer(IK), allocatable          :: res(:)
    allocate (res, source=this%s)
  end function mobbrmsd_state_dump
!
!| dump header as integer array (for python interface api)
  pure function mobbrmsd_state_dump_real(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    real(RK), allocatable             :: res(:)
    allocate (res, source=this%z)
  end function mobbrmsd_state_dump_real
!
!| load integer array as header (for python interface api)
  pure subroutine mobbrmsd_state_load(this, s, z)
    type(mobbrmsd_state), intent(inout) :: this
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

