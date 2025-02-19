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
  public :: mobbrmsd_state_bbgap
  public :: mobbrmsd_state_autovariance
  public :: mobbrmsd_state_squared_deviation
  public :: mobbrmsd_state_mean_squared_deviation
  public :: mobbrmsd_state_rmsd
  public :: mobbrmsd_state_upperbound_as_rmsd
  public :: mobbrmsd_state_lowerbound_as_rmsd
  public :: mobbrmsd_state_n_eval
  public :: mobbrmsd_state_log_n_eval
  public :: mobbrmsd_state_log_sum_n_eval
  public :: mobbrmsd_state_eval_ratio
  public :: mobbrmsd_state_log_eval_ratio
  public :: mobbrmsd_state_has_rotation_matrix
  public :: mobbrmsd_state_is_finished
  public :: mobbrmsd_state_rotation_matrix
  public :: mobbrmsd_state_dump
  public :: mobbrmsd_state_dump_real
  public :: mobbrmsd_state_dump_rotation
  public :: mobbrmsd_state_load
  public :: mobbrmsd_state_destroy
  public :: mobbrmsd_state_attributes
  public :: mobbrmsd_state_RECIPROCAL_OF_N
  public :: mobbrmsd_state_INDEX_TO_AUTOCORR
  public :: mobbrmsd_state_INDEX_TO_UPPERBOUND
  public :: mobbrmsd_state_INDEX_TO_LOWERBOUND
  public :: mobbrmsd_state_INDEX_TO_N_EVAL
  public :: mobbrmsd_state_INDEX_TO_LOG_RATIO
  public :: mobbrmsd_state_INDEX_TO_ROTMAT
  public :: mobbrmsd_state_NOT_YET_RUN_FLAG
  public :: mobbrmsd_state_IS_FINISHED_FLAG
  public :: mobbrmsd_state_ISNOT_FINISHED_FLAG
!&<
  integer(IK), parameter :: mobbrmsd_state_RECIPROCAL_OF_N      = 1
  !! Index to auto correlation
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
!
  integer(IK), parameter :: mobbrmsd_state_NOT_YET_RUN_FLAG     = -1
  !! Flag as not yet running
  integer(IK), parameter :: mobbrmsd_state_IS_FINISHED_FLAG     = 0
  !! Flag as finished
  integer(IK), parameter :: mobbrmsd_state_ISNOT_FINISHED_FLAG  = 1
  !! Flag as incomplete
!&>
!| mobbrmsd_state
  type mobbrmsd_state
    sequence
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
   &   RN => mobbrmsd_state_RECIPROCAL_OF_N, &
   &   AC => mobbrmsd_state_INDEX_TO_AUTOCORR, &
   &   UB => mobbrmsd_state_INDEX_TO_UPPERBOUND, &
   &   LB => mobbrmsd_state_INDEX_TO_LOWERBOUND, &
   &   NE => mobbrmsd_state_INDEX_TO_N_EVAL, &
   &   LR => mobbrmsd_state_INDEX_TO_LOG_RATIO, &
   &   RT => mobbrmsd_state_INDEX_TO_ROTMAT &
    )
      this%s = s
      allocate (z(RT + MAX(d, 0)**2 - 1))
      z(RN) = ONE / n
      z(AC) = ZERO
      z(UB) = RHUGE
      z(LB) = -RHUGE
      z(NE) = -RHUGE
      z(LR) = ZERO
      if (d > 0) call eye(d, z(RT))
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
    if (SIZE(this%z) >= mobbrmsd_state_INDEX_TO_UPPERBOUND) then
      res = this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND)
    else
      res = RHUGE
    end if
  end function mobbrmsd_state_upperbound
!
!| returns lowerbound
  pure elemental function mobbrmsd_state_lowerbound(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    if (SIZE(this%z) >= mobbrmsd_state_INDEX_TO_LOWERBOUND) then
      res = this%z(mobbrmsd_state_INDEX_TO_LOWERBOUND)
    else
      res = -RHUGE
    end if
  end function mobbrmsd_state_lowerbound
!
!| returns upperbound - lowerbound
  pure elemental function mobbrmsd_state_bbgap(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    if (SIZE(this%z) >= mobbrmsd_state_INDEX_TO_LOWERBOUND) then
      res = this%z(mobbrmsd_state_INDEX_TO_UPPERBOUND) &
     &    - this%z(mobbrmsd_state_INDEX_TO_LOWERBOUND)
    else
      res = -RHUGE
    end if
  end function mobbrmsd_state_bbgap
!
!| returns autovariance
  pure elemental function mobbrmsd_state_autovariance(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                         :: res
    if (SIZE(this%z) >= mobbrmsd_state_INDEX_TO_AUTOCORR) then
      res = this%z(mobbrmsd_state_INDEX_TO_AUTOCORR)
    else
      res = ZERO
    end if
  end function mobbrmsd_state_autovariance
!
!| returns squared deviation
  pure elemental function mobbrmsd_state_squared_deviation(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  ac => mobbrmsd_state_INDEX_TO_AUTOCORR, &
   &  ub => mobbrmsd_state_INDEX_TO_UPPERBOUND &
   &  )
    if (SIZE(this%z) >= UB) then
      res = MAX(ZERO, this%z(AC) + TWO * this%z(UB))
    else
      res = RHUGE
    end if
    end associate
  end function mobbrmsd_state_squared_deviation
!
!| returns squared deviation
  pure elemental function mobbrmsd_state_mean_squared_deviation(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  rn => mobbrmsd_state_RECIPROCAL_OF_N, &
   &  ac => mobbrmsd_state_INDEX_TO_AUTOCORR, &
   &  ub => mobbrmsd_state_INDEX_TO_UPPERBOUND &
   &  )
    if (SIZE(this%z) >= UB) then
      res = MAX(ZERO, this%z(RN) * (this%z(AC) + TWO * this%z(UB)))
    else
      res = RHUGE
    end if
    end associate
  end function mobbrmsd_state_mean_squared_deviation
!
!| returns rmsd
  pure elemental function mobbrmsd_state_rmsd(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  rn => mobbrmsd_state_RECIPROCAL_OF_N, &
   &  ac => mobbrmsd_state_INDEX_TO_AUTOCORR, &
   &  ub => mobbrmsd_state_INDEX_TO_UPPERBOUND &
   &  )
    if (SIZE(this%z) >= UB) then
      res = SQRT(MAX(ZERO, this%z(RN) * (this%z(AC) + TWO * this%z(UB))))
    else
      res = RHUGE
    end if
    end associate
  end function mobbrmsd_state_rmsd
!
!| returns lowerbound as rmsd
  pure elemental function mobbrmsd_state_upperbound_as_rmsd(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  rn => mobbrmsd_state_RECIPROCAL_OF_N, &
   &  ac => mobbrmsd_state_INDEX_TO_AUTOCORR, &
   &  ub => mobbrmsd_state_INDEX_TO_UPPERBOUND &
   &  )
    if (SIZE(this%z) >= ub) then
      res = SQRT(MAX(ZERO, this%z(rn) * (this%z(ac) + TWO * this%z(ub))))
    else
      res = ZERO
    end if
    end associate
  end function mobbrmsd_state_upperbound_as_rmsd
!
!| returns lowerbound as rmsd
  pure elemental function mobbrmsd_state_lowerbound_as_rmsd(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                          :: res
    associate (&
   &  rn => mobbrmsd_state_RECIPROCAL_OF_N, &
   &  ac => mobbrmsd_state_INDEX_TO_AUTOCORR, &
   &  lb => mobbrmsd_state_INDEX_TO_LOWERBOUND &
   &  )
    if (SIZE(this%z) >= lb) then
      res = SQRT(MAX(ZERO, this%z(rn) * (this%z(ac) + TWO * this%z(lb))))
    else
      res = ZERO
    end if
    end associate
  end function mobbrmsd_state_lowerbound_as_rmsd
!
!| returns log_n_eval
  pure elemental function mobbrmsd_state_log_n_eval(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                         :: res
    associate (NE => mobbrmsd_state_INDEX_TO_N_EVAL)
    if (SIZE(this%z) >= NE) then
      res = LOG(this%z(NE))
    else
      res = -HUGE(0.0_RK)
    end if
    end associate
  end function mobbrmsd_state_log_n_eval
!
!| returns n_eval
  pure elemental function mobbrmsd_state_n_eval(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    integer(IK)                       :: res
    associate (NE => mobbrmsd_state_INDEX_TO_N_EVAL)
    if (SIZE(this%z) >= NE) then
      res = NINT(this%z(NE), IK)
    else
      res = 0_IK
    end if
    end associate
  end function mobbrmsd_state_n_eval
!
  pure function mobbrmsd_state_log_sum_n_eval(n, state) result(res)
    integer(IK), intent(in)          :: n
    type(mobbrmsd_state), intent(in) :: state(n)
    real(RK)                         :: maxn, res
    integer(IK)                      :: i
    if (n < 1) then
      res = -HUGE(0.0_RK)
      return
    end if
    associate (NE => mobbrmsd_state_INDEX_TO_N_EVAL)
      maxn = state(1)%z(NE)
      do i = 2, n
        maxn = MAX(maxn, state(i)%z(NE))
      end do
      maxn = LOG(maxn)
      res = 0.0_RK
      do i = 1, n
        res = res + EXP(LOG(state(i)%z(NE)) - maxn)
      end do
      res = LOG(res) + maxn
    end associate
  end function mobbrmsd_state_log_sum_n_eval
!
!| returns eval_ratio
  pure elemental function mobbrmsd_state_eval_ratio(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK)                         :: res
    associate (LR => mobbrmsd_state_INDEX_TO_LOG_RATIO)
    if (SIZE(this%z) >= LR) then
      res = EXP(this%z(LR))
    else
      res = ZERO
    end if
    end associate
  end function mobbrmsd_state_eval_ratio
!
!| rotation matrix
  pure elemental function mobbrmsd_state_has_rotation_matrix(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    logical                          :: res
    associate (rt => mobbrmsd_state_INDEX_TO_ROTMAT)
      res = SIZE(this%z) >= rt
    end associate
  end function mobbrmsd_state_has_rotation_matrix
!
!| rotation matrix
  pure subroutine mobbrmsd_state_rotation_matrix(this, R)
    type(mobbrmsd_state), intent(in) :: this
    !! this
    real(RK), intent(inout)          :: R(*)
    !! coordinate
    integer                          :: n
    if (.not. mobbrmsd_state_has_rotation_matrix(this)) return
    associate (rt => mobbrmsd_state_INDEX_TO_ROTMAT)
      n = SIZE(this%z(rt:))
      R(:n) = this%z(rt:)
    end associate
  end subroutine mobbrmsd_state_rotation_matrix
!
!| returns log_eval_ratio
  pure elemental function mobbrmsd_state_log_eval_ratio(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    real(RK)                         :: res
    if (SIZE(this%z) >= mobbrmsd_state_INDEX_TO_LOG_RATIO) then
      res = this%z(mobbrmsd_state_INDEX_TO_LOG_RATIO)
    else
      res = -RHUGE
    end if
  end function mobbrmsd_state_log_eval_ratio
!
!| Returns bb process is finished.
  pure elemental function mobbrmsd_state_is_finished(this) result(res)
    type(mobbrmsd_state), intent(in)  :: this
    !! mobbrmsd_state
    logical                           :: res
    integer                           :: ss
    ss = SIZE(this%s)
    res = this%s(ss) == mobbrmsd_state_IS_FINISHED_FLAG
  end function mobbrmsd_state_is_finished
!
!| dump header as integer array (for python interface api)
  pure function mobbrmsd_state_dump(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    integer(IK), allocatable         :: res(:)
    allocate (res, source=this%s)
  end function mobbrmsd_state_dump
!
!| dump header as integer array (for python interface api)
  pure function mobbrmsd_state_dump_real(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    real(RK), allocatable            :: res(:)
    associate ( &
   &   RN => mobbrmsd_state_RECIPROCAL_OF_N, &
   &   AC => mobbrmsd_state_INDEX_TO_AUTOCORR, &
   &   UB => mobbrmsd_state_INDEX_TO_UPPERBOUND, &
   &   LB => mobbrmsd_state_INDEX_TO_LOWERBOUND, &
   &   NE => mobbrmsd_state_INDEX_TO_N_EVAL, &
   &   LR => mobbrmsd_state_INDEX_TO_LOG_RATIO &
    )
    if (SIZE(this%z) >= LR) then
      allocate (res, source=this%z(:LR))
    else
      allocate (res(LR))
      res(RN) = ONE
      res(AC) = ZERO
      res(UB) = RHUGE
      res(LB) = -RHUGE
      res(NE) = -RHUGE
      res(LR) = ZERO
    end if
    end associate
  end function mobbrmsd_state_dump_real
!
!| dump header as integer array (for python interface api)
  pure function mobbrmsd_state_dump_rotation(this) result(res)
    type(mobbrmsd_state), intent(in) :: this
    !! mobbrmsd_header
    real(RK), allocatable            :: res(:)
    associate ( &
   &   RT => mobbrmsd_state_INDEX_TO_ROTMAT &
    )
    if (mobbrmsd_state_has_rotation_matrix(this)) then
      allocate (res, source=this%z(RT:))
    else
      allocate (res(0))
    end if
    end associate
  end function mobbrmsd_state_dump_rotation
!
!| load integer array as header (for python interface api)
  pure subroutine mobbrmsd_state_load(this, s, z, r)
    type(mobbrmsd_state), intent(inout) :: this
    !! mobbrmsd_header
    integer(IK), intent(in)             :: s(:)
    !! state integer array
    real(RK), intent(in), optional      :: z(:)
    !! state real array
    real(RK), intent(in), optional      :: r(:)
    !! rotation
    integer(IK), allocatable            :: s_(:)
    real(RK), allocatable               :: z_(:)
    allocate (s_, source=s)
    call MOVE_ALLOC(from=s_, to=this%s)
    if (PRESENT(z)) then
      allocate (z_, source=z)
      if (PRESENT(r)) z_ = [z_, r]
    else
      allocate (z_(0))
    end if
    call MOVE_ALLOC(from=z_, to=this%z)
  end subroutine mobbrmsd_state_load
!
!| attributes
  pure subroutine mobbrmsd_state_attributes(d, s, n_int, n_float)
    integer(IK), intent(in)            :: d, s(:)
    integer(IK), intent(out), optional :: n_int
    !! length of state array s
    integer(IK), intent(out), optional :: n_float
    !! length of state array z
    if (PRESENT(n_int)) n_int = SIZE(s)
    if (PRESENT(n_float)) n_float = mobbrmsd_state_INDEX_TO_ROTMAT - 1
  end subroutine mobbrmsd_state_attributes
!
  pure elemental subroutine mobbrmsd_state_destroy(this)
    type(mobbrmsd_state), intent(inout) :: this
    if (ALLOCATED(this%s)) deallocate (this%s)
    if (ALLOCATED(this%z)) deallocate (this%z)
  end subroutine mobbrmsd_state_destroy
end module mod_mobbrmsd_state

