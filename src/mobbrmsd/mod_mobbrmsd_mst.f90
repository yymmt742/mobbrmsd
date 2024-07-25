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
  type prim_list
    sequence
    integer(IK) :: i
    integer(IK) :: j
    integer(IK) :: p
    real(RK)    :: c
  end type prim_list
!
contains
!
! !| nearest_neighbor calculation
! subroutine mobbrmsd_nearest_neighbor( &
!&             n_target, header, state, &
!&             W, &
!&             cutoff, difflim, maxeval, &
!&             remove_com, sort_by_g, &
!&             mask, nnval, nnidx &
!&            )
!   integer(IK), intent(in)             :: n_target
!   !! number of target coordinates
!   type(mobbrmsd), intent(in)          :: header
!   !! mobbrmsd_header
!   type(mobbrmsd_state), intent(inout) :: state(n_target)
!   !! mobbrmsd_state, the result is contained in this structure.
!   real(kind=RK), intent(inout)        :: W(*)
!  !! work memory, must be larger than header%memsize() * mobbrmsd_num_threads()
!   real(RK), intent(in), optional      :: cutoff
!   !! The search ends when lowerbound is determined to be greater than to cutoff.
!   real(RK), intent(in), optional      :: difflim
!   !! The search ends when the difference between the lower and upper bounds is less than difflim.
!   integer(IK), intent(in), optional   :: maxeval
!   !! The search ends when ncount exceeds maxiter.
!   logical, intent(in), optional       :: remove_com
!   !! if true, remove centroids. default [.true.]
!   logical, intent(in), optional       :: sort_by_g
!   !! if true, row is sorted respect to G of reference coordinate. default [.true.]
!   logical, intent(in), optional       :: mask(n_target)
!   !! If .false., skip the calculation.
!   real(RK), intent(out), optional     :: nnval
!   !! nearest neighbor index
!   integer(IK), intent(out), optional  :: nnidx
!   !! nearest neighbor index
!   type(priority_list), allocatable    :: pl(:)
!   real(kind=RK)                       :: cutoff_global, ub
!   integer(kind=IK)                    :: i, j, itgt, ntgt, ypnt, wpnt, ldy, ldw
!
!   ldy = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
!   ldw = mobbrmsd_memsize(header)
!
!   if (PRESENT(cutoff)) then
!     cutoff_global = MERGE(RHUGE, cutoff, cutoff < ZERO)
!   else
!     cutoff_global = RHUGE
!   end if
!
!   if (PRESENT(mask)) then
!     cutoff_global = MIN(MINVAL(mobbrmsd_state_rmsd(state), mask), cutoff_global)
!     ntgt = COUNT(mask)
!   else
!     cutoff_global = MIN(MINVAL(mobbrmsd_state_rmsd(state)), cutoff_global)
!     ntgt = n_target
!   end if
!
!   allocate (pl(ntgt))
!
!   if (PRESENT(mask)) then
!     j = 0
!     do i = 1, n_target
!       if (.not. mask(i)) cycle
!       j = j + 1
!       pl(j) = priority_list(i, mobbrmsd_state_upperbound(state(i)))
!     end do
!   else
!     do concurrent(i=1:ntgt)
!       pl(i) = priority_list(i, mobbrmsd_state_upperbound(state(i)))
!     end do
!   end if
!
!   if (ntgt > 1) call qs(ntgt, pl)
!
!   i = 0
!
!   !$omp parallel private(itgt, ub)
!   do
!     !$omp critical
!     i = i + 1
!     itgt = i
!     ub = cutoff_global
!     !$omp end critical
!
!     if (itgt > ntgt) exit
!
!     itgt = pl(itgt)%i
!
!     if (mobbrmsd_is_finished(header, state(itgt))) cycle
!     if (ub < mobbrmsd_state_lowerbound_as_rmsd(state(itgt))) cycle
!
!     wpnt = (itgt - 1) * ldw + 1
!
!     call mobbrmsd_restart( &
!    &       header, &
!    &       state(itgt), &
!    &       W(wpnt), &
!    &       cutoff=ub, &
!    &       difflim=difflim, &
!    &       maxeval=maxeval &
!    &      )
!     ub = mobbrmsd_state_rmsd(state(itgt))
!
!     !$omp critical
!     cutoff_global = MIN(ub, cutoff_global)
!     !$omp end critical
!   end do
!   !$omp end parallel
!
!   if (PRESENT(mask)) then
!     if (PRESENT(nnval)) nnval = MINVAL(mobbrmsd_state_rmsd(state), mask)
!     if (PRESENT(nnidx)) nnidx = MINLOC(mobbrmsd_state_rmsd(state), 1, mask)
!   else
!     if (PRESENT(nnval)) nnval = MINVAL(mobbrmsd_state_rmsd(state))
!     if (PRESENT(nnidx)) nnidx = MINLOC(mobbrmsd_state_rmsd(state), 1)
!   end if
!
! end subroutine mobbrmsd_nearest_neighbor
!
!| min_span_tree construction
! subroutine mobbrmsd_min_span_tree( &
!&             n_target, &
!&             header, &
!&             state, &
!&             W, &
!&             cutoff, &
!&             difflim, &
!&             maxeval, &
!&             remove_com, &
!&             sort_by_g, &
!&             edges, &
!&             weights &
!&            )
!   type(mobbrmsd), intent(in)          :: header
!   !! mobbrmsd_header
!   integer(IK), intent(in)             :: n_target
!   !! number of target coordinates
!   type(mobbrmsd_state), intent(inout) :: state(*)
!   !! mobbrmsd_state, the result is contained in this structure.
!   real(kind=RK), intent(inout)        :: W(*)
!   !! work memory, must be larger than header%memsize() * n_target * (n_target-1) / 2
!   real(RK), intent(in), optional      :: cutoff
!   !! The search ends when lowerbound is determined to be greater than to cutoff.
!   real(RK), intent(in), optional      :: difflim
!   !! The search ends when the difference between the lower and upper bounds is less than difflim.
!   integer(IK), intent(in), optional   :: maxeval
!   !! The search ends when ncount exceeds maxiter.
!   logical, intent(in), optional       :: remove_com
!   !! if true, remove centroids. default [.true.]
!   logical, intent(in), optional       :: sort_by_g
!   !! if true, row is sorted respect to G of reference coordinate. default [.true.]
!   integer(IK), intent(inout)          :: edges(2, n_target)
!   !! mst edges
!   real(RK), intent(inout)             :: weights(n_target)
!   !! mst weights
!   type(prim_list)                     :: pl(n_target)
!   real(kind=RK)                       :: cutoff_global, ub
!   integer(kind=IK)                    :: i, j, itgt, ntgt, ypnt, wpnt, ldy, ldw
!
!   ldy = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
!   ldw = mobbrmsd_memsize(header)
!
!   if (PRESENT(cutoff)) then
!     cutoff_global = MERGE(RHUGE, cutoff, cutoff < ZERO)
!   else
!     cutoff_global = RHUGE
!   end if
!
!   if (PRESENT(mask)) then
!     cutoff_global = MIN(MINVAL(mobbrmsd_state_rmsd(state), mask), cutoff_global)
!     ntgt = COUNT(mask)
!   else
!     cutoff_global = MIN(MINVAL(mobbrmsd_state_rmsd(state)), cutoff_global)
!     ntgt = n_target
!   end if
!
!   if (PRESENT(mask)) then
!     j = 0
!     do i = 1, n_target
!       if (.not. mask(i)) cycle
!       j = j + 1
!       pl(j) = priority_list(i, mobbrmsd_state_upperbound(state(i)))
!     end do
!   else
!     do concurrent(i=1:ntgt)
!       pl(i) = priority_list(i, mobbrmsd_state_upperbound(state(i)))
!     end do
!   end if
!
!   if (ntgt > 1) call qs(ntgt, pl)
!
!   i = 0
!
!   !$omp parallel private(itgt, ub)
!   do
!     !$omp critical
!     i = i + 1
!     itgt = i
!     ub = cutoff_global
!     !$omp end critical
!
!     if (itgt > ntgt) exit
!
!     itgt = pl(itgt)%i
!
!     if (mobbrmsd_is_finished(header, state(itgt))) cycle
!     if (ub < mobbrmsd_state_lowerbound_as_rmsd(state(itgt))) cycle
!
!     wpnt = (itgt - 1) * ldw + 1
!
!     call mobbrmsd_restart( &
!    &  header, &
!    &  state(itgt), &
!    &  W(wpnt), &
!    &  cutoff=ub, &
!    &  difflim=difflim, &
!    &  maxeval=maxeval, &
!    &      )
!     ub = mobbrmsd_state_rmsd(state(itgt))
!
!     !$omp critical
!     cutoff_global = MIN(ub, cutoff_global)
!     !$omp end critical
!   end do
!   !$omp end parallel
!
!   if (PRESENT(mask)) then
!     if (PRESENT(nnval)) nnval = MINVAL(mobbrmsd_state_rmsd(state), mask)
!     if (PRESENT(nnidx)) nnidx = MINLOC(mobbrmsd_state_rmsd(state), 1, mask)
!   else
!     if (PRESENT(nnval)) nnval = MINVAL(mobbrmsd_state_rmsd(state))
!     if (PRESENT(nnidx)) nnidx = MINLOC(mobbrmsd_state_rmsd(state), 1)
!   end if
!
! end subroutine mobbrmsd_nearest_neighbor
!
!| minimum spanning tree construction.
  subroutine mobbrmsd_min_span_tree( &
 &             n_target, &
 &             header, &
 &             state, &
 &             X, &
 &             W, &
 &             cutoff, &
 &             difflim, &
 &             maxeval, &
 &             remove_com, &
 &             sort_by_g, &
 &             edges, &
 &             weights &
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
!   !! work memory, must be larger than header%memsize() * n_target * (n_target-1) / 2
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
    type(prim_list)                     :: pl(n_target)
    real(RK)                            :: rmsd, ub, cutoff_global
    integer(kind=IK)                    :: i, j, k, itgt, xpnt, ypnt, wpnt, l0, ldx, ldw
!
    l0 = n_target * (n_target - 1) / 2
    ldx = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
    ldw = mobbrmsd_memsize(header)
!
    !$omp parallel do private(itgt, xpnt, ypnt, wpnt)
    do j = 1, n_target
      do i = 1, j - 1
        call cantor_pair(i, j, itgt)
        xpnt = (i - 1) * ldx + 1
        ypnt = (j - 1) * ldx + 1
        wpnt = (itgt - 1) * ldw + 1
        call mobbrmsd_run( &
       &       header, &
       &       state(itgt), &
       &       X(xpnt), &
       &       X(ypnt), &
       &       W(wpnt), &
       &       cutoff=cutoff, &
       &       difflim=difflim, &
       &       maxeval=1, &
       &       remove_com=remove_com, &
       &       sort_by_g=sort_by_g  &
       &      )
      end do
    end do
    !$omp end parallel do
!
    if (PRESENT(cutoff)) then
      cutoff_global = MERGE(RHUGE, cutoff, cutoff < ZERO)
    else
      cutoff_global = RHUGE
    end if
!
    do k = 1, l0
      call cantor_pair_inverse(k, i, j)
    end do
    do concurrent(i=1:n_target)
      pl(i) = prim_list(i, 0, 0, RHUGE)
    end do
    k = MINLOC(mobbrmsd_state_lowerbound_as_rmsd(state(:l0)), 1)
    call cantor_pair_inverse(k, i, j)
    call prim_list_swap(pl(1), pl(i))
!
    do k = 2, n_target
      do concurrent(i=1:k - 1)
        pl(i)%c = RHUGE
      end do
      ub = cutoff_global
!
      do i = 1, k - 1
        do j = k, n_target
          call cantor_pair(pl(i)%i, pl(j)%i, itgt)
          if (ub < mobbrmsd_state_lowerbound_as_rmsd(state(itgt))) cycle
          wpnt = (itgt - 1) * ldw + 1
          if (.not. mobbrmsd_is_finished(header, state(itgt))) then
            call mobbrmsd_restart( &
           &       header, &
           &       state(itgt), &
           &       W(wpnt), &
           &       cutoff=ub, &
           &       difflim=difflim, &
           &       maxeval=maxeval &
           &      )
          end if
          rmsd = mobbrmsd_state_rmsd(state(itgt))
          ub = MIN(ub, rmsd)
          if (rmsd < pl(i)%c) then
            pl(i)%j = pl(j)%i
            pl(i)%p = j
            pl(i)%c = rmsd
          end if
        end do
      end do
      j = MINLOC(pl(:k - 1)%c, 1)
      if (PRESENT(edges)) then
        edges(1, k - 1) = pl(j)%i
        edges(2, k - 1) = pl(j)%j
      end if
      if (PRESENT(weights)) then
        weights(k - 1) = pl(j)%c
      end if
      call prim_list_swap(pl(k), pl(pl(j)%p))
    end do
!
  end subroutine mobbrmsd_min_span_tree
!
! ---
!
! pure recursive subroutine qs(n, q)
!   integer(IK), intent(in)        :: n
!   type(prim_list), intent(inout) :: q(*)
!   type(prim_list)                :: t
!   integer(IK)                    :: i, j, p
!   j = n; i = 1; p = j / 2
!   do
!     do while (q(i)%u < q(p)%u); i = i + 1; end do
!     do while (q(j)%u > q(p)%u); j = j - 1; end do
!     if (i >= j) exit
!     t = q(i); q(i) = q(j); q(j) = t
!     if (i == p) then; p = j
!     elseif (j == p) then; p = i
!     end if
!     i = i + 1; j = j - 1
!   end do
!   if (2 < i) call qs(i - 1, q(1))
!   if (j < n - 1) call qs(n - j, q(j + 1))
! end subroutine qs
!
  pure elemental subroutine prim_list_swap(a, b)
    type(prim_list), intent(inout) :: a, b
    type(prim_list)                :: t
    t = b
    b = a
    a = t
  end subroutine prim_list_swap
!
  pure elemental subroutine cantor_pair(i, j, k)
    integer(IK), intent(in)    :: i, j
    integer(IK), intent(inout) :: k
    if (i > j) then
      k = (i - 2) * (i - 1) / 2 + j
    else
      k = (j - 2) * (j - 1) / 2 + i
    end if
  end subroutine cantor_pair
!
  pure elemental subroutine cantor_pair_inverse(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    j = INT((SQRT(real(8 * k - 7, RK)) - ONE), IK) / 2 + 2 ! cantor_pair_inverse
    i = k - (j - 1) * (j - 2) / 2
  end subroutine cantor_pair_inverse
!
end module mod_mobbrmsd_mst

