!| Configure a minimum global tree with mobbrmsd. <br>
!  MST construction is based on the Prim's algorithm,
!  but the useless calculations are reduced using the cutoff possibilities of mobbrmsd.
module mod_mobbrmsd_mst
!$ use omp_lib
  use mod_params, only: &
 &      IK, &
 &      RK, &
 &      ONE => RONE, &
 &      ZERO => RZERO, &
 &      RHUGE
  use mod_mobbrmsd_header
  use mod_mobbrmsd_state
  use mod_mobbrmsd
  use mod_forbar
  use mod_forbar_collections
  implicit none
  private
  public :: mobbrmsd_min_span_tree
!
  type priority_list
    sequence
    integer(IK) :: i
    real(RK)    :: u
  end type priority_list
!
contains
!
  !| nearest_neighbor calculation
  subroutine mobbrmsd_nearest_neighbor(n_target, header, state, X, Y, W, &
 &                                     cutoff, difflim, maxeval, &
 &                                     mask, nnval, nnidx)
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
    type(priority_list), allocatable    :: pl(:)
    real(kind=RK)                       :: cutoff_global, ub
    integer(kind=IK)                    :: i, j, itgt, ntgt, ypnt, wpnt, ldy, ldw
!
    ldy = header%n_dims() * header%n_atoms()
    ldw = header%memsize()
!
    if (PRESENT(cutoff)) then
      cutoff_global = MERGE(RHUGE, cutoff, cutoff < ZERO)
    else
      cutoff_global = RHUGE
    end if
!
    if (PRESENT(mask)) then
      cutoff_global = MIN(MINVAL(state%upperbound(), mask), cutoff_global)
      ntgt = COUNT(mask)
    else
      cutoff_global = MIN(MINVAL(state%upperbound()), cutoff_global)
      ntgt = n_target
    end if
!
    allocate (pl(ntgt))
!
    if (PRESENT(mask)) then
      j = 0
      do i = 1, n_target
        if (.not. mask(i)) cycle
        j = j + 1
        pl(j) = priority_list(i, state(i)%upperbound())
      end do
    else
      do concurrent(i=1:ntgt)
        pl(i) = priority_list(i, state(i)%upperbound())
      end do
    end if
!
    if (ntgt > 1) call qs(ntgt, pl)
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
      if (itgt > ntgt) exit
!
      itgt = pl(itgt)%i
!
      if (bb_list_is_finished(header%q, state(itgt)%s)) cycle
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
      ub = state(itgt)%rmsd()
!
      !$omp critical
      cutoff_global = MIN(ub, cutoff_global)
      !$omp end critical
    end do
    !$omp end parallel
!
    if (PRESENT(mask)) then
      if (PRESENT(nnval)) nnval = MINVAL(state%upperbound(), mask)
      if (PRESENT(nnidx)) nnidx = MINLOC(state%upperbound(), 1, mask)
    else
      if (PRESENT(nnval)) nnval = MINVAL(state%upperbound())
      if (PRESENT(nnidx)) nnidx = MINLOC(state%upperbound(), 1)
    end if
!
  end subroutine mobbrmsd_nearest_neighbor
!
!| min_span_tree construction
  subroutine mobbrmsd_min_span_tree(n_target, header, state, X, W, &
 &                                  cutoff, difflim, maxeval, &
 &                                  edges, weights, show_progress, &
 &                                  verbose &
 )
    integer(IK), intent(in)             :: n_target
    !! number of coordinates
    type(mobbrmsd_header), intent(in)   :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout) :: state(n_target, n_target)
    !! mobbrmsd_state, the result is contained in this structure.
    real(kind=RK), intent(in)           :: X(*)
    !! coordinate sequence
    real(kind=RK), intent(inout)        :: W(*)
    !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
    real(RK), intent(in), optional      :: cutoff
    !! The search ends when lowerbound is determined to be greater than to cutoff.
    real(RK), intent(in), optional      :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional   :: maxeval
    !! The search ends when ncount exceeds maxiter.
    integer(IK), intent(out), optional  :: edges(2, n_target - 1)
    !! minimum spanning tree edges
    real(RK), intent(out), optional     :: weights(n_target - 1)
    !! minimum spanning tree weights
    logical, intent(in), optional       :: show_progress
    !! if true, show progress bar
    logical, intent(in), optional       :: verbose
    !! show progress bar
    type(forbar)                        :: fbar
    logical                             :: verbose_
    logical                             :: mask(n_target)
    integer(kind=IK)                    :: list(2, n_target)
    real(RK)                            :: vval(n_target - 1), nnval, cutoff_
    integer(kind=IK)                    :: i, j, xpnt, ldx, nnidx
!
    ldx = header%n_dims() * header%n_atoms()
!
!   Initialize header
    do concurrent(i=1:n_target, j=1:n_target)
      state(i, j) = mobbrmsd_state(header)
    end do
!
    mask(:) = .true.
    mask(1) = .false.
    list(1, 1) = 0
    list(2, 1) = 1
!
    if (PRESENT(verbose)) then
      verbose_ = verbose
    else
      verbose_ = .false.
    end if
!
    if (verbose_) then
      fbar = progress_bar(limit=n_target - 1)
    end if
!
    do j = 1, n_target - 1
!
      if (verbose_) call fbar%update_and_display()
!
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
    if (verbose_) call fbar%clean_up()
!
    do j = 1, n_target
      do i = 1, j - 1
        if (state(i, j)%upperbound() < state(j, i)%upperbound()) then
          state(j, i) = state(i, j)
        else
          state(i, j) = state(j, i)
        end if
      end do
    end do
!
    if (PRESENT(edges)) then
      do concurrent(i=1:n_target - 1)
        edges(1, i) = list(1, i + 1)
        edges(2, i) = list(2, i + 1)
      end do
    end if
!
    if (PRESENT(weights)) then
      do concurrent(i=1:n_target - 1)
        weights(i) = vval(i)
      end do
    end if
!
  end subroutine mobbrmsd_min_span_tree
!
! ---
!
  pure recursive subroutine qs(n, q)
    integer, intent(in)                :: n
    type(priority_list), intent(inout) :: q(*)
    type(priority_list)                :: t
    integer                            :: i, j, p
    j = n; i = 1; p = j / 2
    do
      do while (q(i)%u < q(p)%u); i = i + 1; end do
      do while (q(j)%u > q(p)%u); j = j - 1; end do
      if (i >= j) exit
      t = q(i); q(i) = q(j); q(j) = t
      if (i == p) then; p = j
      elseif (j == p) then; p = i
      end if
      i = i + 1; j = j - 1
    end do
    if (2 < i) call qs(i - 1, q(1))
    if (j < n - 1) call qs(n - j, q(j + 1))
  end subroutine qs
!
end module mod_mobbrmsd_mst

