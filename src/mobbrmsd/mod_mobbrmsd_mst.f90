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
  use mod_mobbrmsd_state
  use mod_mobbrmsd
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
  subroutine mobbrmsd_nearest_neighbor( &
 &             n_target, header, state, &
 &             W, &
 &             cutoff, difflim, maxeval, &
 &             remove_com, sort_by_g, &
 &             mask, nnval, nnidx &
 &            )
    integer(IK), intent(in)             :: n_target
    !! number of target coordinates
    type(mobbrmsd), intent(in)          :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout) :: state(n_target)
    !! mobbrmsd_state, the result is contained in this structure.
    real(kind=RK), intent(inout)        :: W(*)
   !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
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
    ldy = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
    ldw = mobbrmsd_memsize(header)
!
    if (PRESENT(cutoff)) then
      cutoff_global = MERGE(RHUGE, cutoff, cutoff < ZERO)
    else
      cutoff_global = RHUGE
    end if
!
    if (PRESENT(mask)) then
      cutoff_global = MIN(MINVAL(mobbrmsd_state_rmsd(state), mask), cutoff_global)
      ntgt = COUNT(mask)
    else
      cutoff_global = MIN(MINVAL(mobbrmsd_state_rmsd(state)), cutoff_global)
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
        pl(j) = priority_list(i, mobbrmsd_state_upperbound(state(i)))
      end do
    else
      do concurrent(i=1:ntgt)
        pl(i) = priority_list(i, mobbrmsd_state_upperbound(state(i)))
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
      if (mobbrmsd_is_finished(header, state(itgt))) cycle
      if (ub < mobbrmsd_state_lowerbound_as_rmsd(state(itgt))) cycle
!
      wpnt = (itgt - 1) * ldw + 1
!
      call mobbrmsd_restart( &
     &  header, &
     &  state(itgt), &
     &  W(wpnt), &
     &  cutoff=ub, &
     &  difflim=difflim, &
     &  maxeval=maxeval, &
     &      )
      ub = mobbrmsd_state_rmsd(state(itgt))
!
      !$omp critical
      cutoff_global = MIN(ub, cutoff_global)
      !$omp end critical
    end do
    !$omp end parallel
!
    if (PRESENT(mask)) then
      if (PRESENT(nnval)) nnval = MINVAL(mobbrmsd_state_rmsd(state), mask)
      if (PRESENT(nnidx)) nnidx = MINLOC(mobbrmsd_state_rmsd(state), 1, mask)
    else
      if (PRESENT(nnval)) nnval = MINVAL(mobbrmsd_state_rmsd(state))
      if (PRESENT(nnidx)) nnidx = MINLOC(mobbrmsd_state_rmsd(state), 1)
    end if
!
  end subroutine mobbrmsd_nearest_neighbor
!
!| min_span_tree construction
  subroutine mobbrmsd_min_span_tree( &
 &             n_target, header, state, &
 &             X, W, &
 &             cutoff, difflim, maxeval, &
 &             remove_com, sort_by_g, &
 &             edges, weights &
 &            )
    integer(IK), intent(in)             :: n_target
    !! number of coordinates
    type(mobbrmsd), intent(in)          :: header
    !! mobbrmsd_header
    type(mobbrmsd_state), intent(inout) :: state(n_target * (n_target - 1) / 2)
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
    logical, intent(in), optional       :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional       :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
    integer(IK), intent(out), optional  :: edges(2, n_target - 1)
    !! minimum spanning tree edges
    real(RK), intent(out), optional     :: weights(n_target - 1)
    !! minimum spanning tree weights
    logical                             :: mask(n_target)
    integer(kind=IK)                    :: list(2, n_target)
    real(RK)                            :: vval(n_target - 1), nnval, cutoff_
    integer(kind=IK)                    :: i, j, itgt, xpnt, ypnt, wpnt, ldx, ldw, nnidx
!
    ldx = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
    ldw = mobbrmsd_memsize(header)
!
    do j = 1, n_target
      do i = 1, j - 1
        !do concurrent(j=1:n_target)
        !  do concurrent(i=1:j - 1)
        itgt = (j - 2) * (j - 1) / 2 + i
        xpnt = (i - 1) * ldx + 1
        ypnt = (j - 1) * ldx + 1
        wpnt = (itgt - 1) * ldw + 1
        call mobbrmsd_run( &
       &  header, state(itgt), &
       &  X(xpnt), X(ypnt), W(wpnt), &
       &  cutoff=cutoff, &
       &  difflim=difflim, &
       &  maxeval=1, &
       &  remove_com=remove_com, &
       &  sort_by_g=sort_by_g  &
       &      )
      end do
    end do
!
    mask(:) = .true.
    mask(1) = .false.
    list(1, 1) = 0
    list(2, 1) = 1
!
    do j = 1, n_target - 1
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
!       call mobbrmsd_nearest_neighbor( &
!      &  n_target, &
!      &  header, &
!      &  state(:, list(2, i)), &
!      &  W, &
!      &  cutoff=cutoff_, &
!      &  difflim=difflim, &
!      &  maxeval=maxeval, &
!      &  remove_com=remove_com, &
!      &  sort_by_g=sort_by_g, &
!      &  mask=mask, &
!      &  nnval=nnval, &
!      &  nnidx=nnidx &
!      & )
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
!   do j = 1, n_target
!     do i = 1, j - 1
!       if (mobbrmsd_state_upperbound(state(i, j)) < mobbrmsd_state_upperbound(state(j, i))) then
!         state(j, i) = state(i, j)
!       else
!         state(i, j) = state(j, i)
!       end if
!     end do
!   end do
!
!   if (PRESENT(edges)) then
!     do concurrent(i=1:n_target - 1)
!       edges(1, i) = list(1, i + 1)
!       edges(2, i) = list(2, i + 1)
!     end do
!   end if
!
!   if (PRESENT(weights)) then
!     do concurrent(i=1:n_target - 1)
!       weights(i) = vval(i)
!     end do
!   end if
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

