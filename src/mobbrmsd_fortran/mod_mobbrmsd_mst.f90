!| Configure a minimum global tree with mobbrmsd. <br>
!  MST construction is based on the Kruscal algorithm,
!  but the useless calculations are reduced using the cutoff possibilities of mobbrmsd.
module mod_mobbrmsd_mst
!$ use omp_lib
  use mod_params, only: &
 &      I8, &
 &      IK, &
 &      RK, &
 &      ONE => RONE, &
 &      ZERO => RZERO, &
 &      RHUGE
  use mod_mobbrmsd_state
  use mod_mobbrmsd
  use mod_mobbrmsd_batch_run
  implicit none
  private
  public :: mobbrmsd_min_span_tree
!
  integer(IK), parameter :: IS_FINISHED = 0
  integer(IK), parameter :: NOT_FINISHED = 1
  integer(IK), parameter :: NOT_EVALUATED = 2
!
  type heap_container
    sequence
    integer(IK) :: p, s
    real(RK)    :: lb, ub
  end type heap_container
!
contains
!| minimum spanning tree construction.
  subroutine mobbrmsd_min_span_tree(n_target&
                                 &, header &
                                 &, X &
                                 &, n_work &
                                 &, remove_com &
                                 &, sort_by_g &
                                 &, verbose &
                                 &, edges &
                                 &, weights &
                                 &, log_n_eval &
                                 &)
    integer(IK), intent(in)             :: n_target
    !! number of coordinates
    type(mobbrmsd), intent(in)          :: header
    !! mobbrmsd_header
    real(RK), intent(in)                :: X(*)
    !! coordinate sequence
    integer(IK), intent(in), optional   :: n_work
    !! Maximum memory storage area size
    logical, intent(in), optional       :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional       :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
    logical, intent(in), optional       :: verbose
    !! if true, enable verbose output
    integer(IK), intent(out), optional  :: edges(2, n_target - 1)
    !! minimum spanning tree edges
    real(RK), intent(out), optional     :: weights(n_target - 1)
    !! minimum spanning tree weights
    real(RK), intent(out), optional     :: log_n_eval
    !! log_n_eval
    type(heap_container)                :: core(n_target - 1) ! heap container
    integer(IK)                         :: par(n_target) ! for union find
    integer(IK)                         :: n_tri, ldw, ldx
    integer(IK)                         :: n_core, n_reco, n_thread, n_edge, n_heap
    integer(IK)                         :: i, j, k
!
    if (.not. PRESENT(edges) .and. .not. PRESENT(weights) .and. .not. PRESENT(log_n_eval)) return
!
    n_edge = n_target - 1
    n_tri = n_target * (n_target - 1) / 2
    ldw = mobbrmsd_memsize(header)
    ldx = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
!
    !$omp parallel
    if (omp_get_thread_num() == 0) n_thread = omp_get_num_threads()
    !$omp end parallel
!
    n_heap = n_tri
    if (PRESENT(n_work)) then
      if (n_work > 0) n_heap = MIN(MAX(n_target * 2, n_work), n_heap)
    end if
    n_heap = MAX(n_heap, n_thread)
!
    block
      type(mobbrmsd_state), allocatable :: state(:)
      integer(IK)          :: hptr(n_heap) ! heap pointer
      integer(IK)          :: sptr(n_heap) ! sorted pointer
      type(heap_container) :: hcon(n_heap) ! heap container
      real(RK)             :: wrkmem(ldw, n_heap)
      real(RK)             :: lambda
!
      allocate (state(n_tri))
!
      call mobbrmsd_batch_tri_run(n_target &
                               &, header &
                               &, state &
                               &, X &
                               &, wrkmem &
                               &, maxeval=0 &
                               &, remove_com=remove_com &
                               &, sort_by_g=sort_by_g &
                               & )
!
      n_core = 0
      n_reco = 0
      lambda = 0.0_RK
      call init_arange(n_target, par)
!
      do while (n_core < n_edge)
!
        call reflesh_heap(n_tri &
                       &, n_target &
                       &, n_edge &
                       &, n_heap &
                       &, n_core &
                       &, lambda &
                       &, par &
                       &, core &
                       &, state &
                       &, n_reco &
                       &, hptr &
                       &, hcon &
                       &)
!
        lambda = hcon(hptr(1))%lb
        !print *, mobbrmsd_state_log_sum_n_eval(n_tri, state)
        do while (n_core < n_edge .and. n_edge - n_core <= n_reco)
          !print *, 'sort_heap'
          call sort_heap(n_edge &
                      &, n_heap &
                      &, n_reco &
                      &, n_core &
                      &, hptr &
                      &, hcon &
                      &, sptr &
                      & )
!
          !print *, 'kruscal'
          call kruscal(n_target &
                    &, n_edge &
                    &, n_heap &
                    &, lambda &
                    &, hcon &
                    &, n_reco &
                    &, n_core &
                    &, par &
                    &, sptr &
                    &, core &
                    & )
          if (n_core == n_edge .or. (n_edge - n_core > n_reco)) exit
!
          !print *, 'update_bounds'
          call update_bounds(ldx &
                          &, ldw &
                          &, n_edge &
                          &, n_core &
                          &, n_heap &
                          &, n_reco &
                          &, header &
                          &, X &
                          &, wrkmem &
                          &, sptr &
                          &, hcon &
                          &, state &
                          &, remove_com &
                          &, sort_by_g &
                          &)
!
          !print *, 'update_heap'
          call update_heap(n_reco, n_heap, hcon, sptr, hptr)
        end do
      end do
!
      if (PRESENT(edges)) then
        do concurrent(k=1:n_edge)
          call cantor_pair_inverse(core(k)%p, i, j)
          edges(1, k) = i
          edges(2, k) = j
        end do
      end if
      if (PRESENT(weights)) then
        do concurrent(k=1:n_edge)
          weights(k) = core(k)%lb
        end do
      end if
      if (PRESENT(log_n_eval)) then
        log_n_eval = mobbrmsd_state_log_sum_n_eval(n_tri, state)
      end if
    end block
!
  end subroutine mobbrmsd_min_span_tree
!
  pure subroutine reflesh_heap(n_tri &
                            &, n_target &
                            &, n_edge &
                            &, n_heap &
                            &, n_core &
                            &, lambda &
                            &, par &
                            &, core &
                            &, state &
                            &, n_reco &
                            &, hptr &
                            &, hcon &
                            &)
    integer(IK), intent(in)             :: n_tri
    integer(IK), intent(in)             :: n_target
    integer(IK), intent(in)             :: n_edge
    integer(IK), intent(in)             :: n_heap
    integer(IK), intent(in)             :: n_core
    real(RK), intent(in)                :: lambda
    integer(IK), intent(in)             :: par(n_edge)
    type(heap_container), intent(in)    :: core(n_core)
    type(mobbrmsd_state), intent(in)    :: state(n_tri)
    integer(IK), intent(inout)          :: n_reco
    integer(IK), intent(inout)          :: hptr(n_heap)
    type(heap_container), intent(inout) :: hcon(n_heap)
    real(RK)                            :: lb
    integer(IK)                         :: i, n_rem
    n_rem = MIN(n_heap, n_tri - n_core)
    i = 0
    n_reco = 0
    call init_arange(n_heap, hptr)
    do while (i < n_tri .and. n_reco < n_heap)
      i = i + 1
      if (is_same(n_target, i, par)) cycle
      n_reco = n_reco + 1
      hcon(hptr(n_reco))%p = i
      hcon(hptr(n_reco))%s = MERGE(IS_FINISHED, NOT_EVALUATED, mobbrmsd_state_is_finished(state(i)))
      hcon(hptr(n_reco))%lb = mobbrmsd_state_lowerbound_as_rmsd(state(i))
      hcon(hptr(n_reco))%ub = mobbrmsd_state_upperbound_as_rmsd(state(i))
    end do
    call make_heap(n_reco, hcon, hptr)
    do while (i < n_tri)
      i = i + 1
      if (is_same(n_target, i, par)) cycle
      lb = mobbrmsd_state_lowerbound_as_rmsd(state(i))
      if (lb < hcon(hptr(1))%lb) then
        hcon(hptr(1))%p = i
        hcon(hptr(1))%s = MERGE(IS_FINISHED, NOT_EVALUATED, mobbrmsd_state_is_finished(state(i)))
        hcon(hptr(1))%lb = lb
        hcon(hptr(1))%ub = mobbrmsd_state_lowerbound_as_rmsd(state(i))
        call down_heap(n_heap, 1, hcon, hptr)
      end if
    end do
  end subroutine reflesh_heap
!
  pure subroutine sort_heap(n_edge &
                         &, n_heap &
                         &, n_reco &
                         &, n_core &
                         &, hptr &
                         &, hcon &
                         &, sptr &
                         & )
    integer(IK), intent(in)          :: n_edge
    integer(IK), intent(in)          :: n_heap
    integer(IK), intent(in)          :: n_reco
    integer(IK), intent(in)          :: n_core
    integer(IK), intent(in)          :: hptr(n_reco)
    type(heap_container), intent(in) :: hcon(n_heap)
    integer(IK), intent(inout)       :: sptr(n_reco)
    integer(IK)                      :: k
    do concurrent(k=1:n_reco)
      sptr(k) = hptr(k)
    end do
    do k = n_reco, 2, -1
      call swap(sptr(1), sptr(k))
      call down_heap(k - 1, 1, hcon, sptr)
    end do
  end subroutine sort_heap
!
  subroutine update_bounds(ldx &
                        &, ldw &
                        &, n_edge &
                        &, n_core &
                        &, n_heap &
                        &, n_reco &
                        &, header &
                        &, X &
                        &, wrkmem &
                        &, sptr &
                        &, hcon &
                        &, state &
                        &, remove_com &
                        &, sort_by_g &
                        & )
    integer(IK), intent(in)             :: ldx
    integer(IK), intent(in)             :: ldw
    integer(IK), intent(in)             :: n_edge
    integer(IK), intent(in)             :: n_core
    integer(IK), intent(in)             :: n_heap
    integer(IK), intent(in)             :: n_reco
    type(mobbrmsd), intent(in)          :: header
    real(RK), intent(in)                :: X(*)
    real(RK), intent(inout)             :: wrkmem(ldw, n_heap)
    integer(IK), intent(in)             :: sptr(n_reco)
    type(heap_container), intent(inout) :: hcon(n_heap)
    type(mobbrmsd_state), intent(inout) :: state(*)
    logical, intent(in), optional       :: remove_com, sort_by_g
    real(RK)                            :: cutoff
    integer(IK)                         :: k, l, n_eval
    !n_eval = MIN(2 * n_edge, n_reco)
    n_eval = MIN(n_edge - n_core, n_reco)
    cutoff = hcon(sptr(n_edge - n_core))%ub * 1.00001_RK
    l = 0
    !$omp parallel private(k)
    do
      !$omp critical
      l = l + 1
      k = l
      !$omp end critical
      if (k > n_eval) exit
      associate (p => hcon(sptr(k))%p&
              &, s => hcon(sptr(k))%s&
              &, lb => hcon(sptr(k))%lb&
              &, ub => hcon(sptr(k))%ub&
              &, st => state(hcon(sptr(k))%p)&
              &, h => sptr(k)&
              & )
        if (s == IS_FINISHED) cycle
        if (s == NOT_EVALUATED) then
          block
            integer(IK) :: i, j, xpnt, ypnt
            call cantor_pair_inverse(p, i, j)
            xpnt = (j - 1) * ldx + 1
            ypnt = (i - 1) * ldx + 1
            call mobbrmsd_run(header &
                           &, st &
                           &, X(xpnt) &
                           &, X(ypnt) &
                           &, wrkmem(1, h) &
                           &, remove_com=remove_com &
                           &, sort_by_g=sort_by_g &
                           & )
          end block
        elseif (s == NOT_FINISHED) then
          call mobbrmsd_restart(header &
                             &, st &
                             &, wrkmem(1, h) &
                             & )
        end if
        lb = mobbrmsd_state_lowerbound_as_rmsd(st)
        ub = mobbrmsd_state_lowerbound_as_rmsd(st)
        s = MERGE(IS_FINISHED, NOT_FINISHED, mobbrmsd_state_is_finished(st))
      end associate
    end do
    !$omp end parallel
  end subroutine update_bounds
!
  pure subroutine kruscal(n_target &
                       &, n_edge &
                       &, n_heap &
                       &, lambda &
                       &, hcon &
                       &, n_reco &
                       &, n_core &
                       &, par &
                       &, sptr &
                       &, core &
                       &)
    integer(IK), intent(in)             :: n_target
    integer(IK), intent(in)             :: n_edge
    integer(IK), intent(in)             :: n_heap
    real(RK), intent(in)                :: lambda
    type(heap_container), intent(in)    :: hcon(n_heap)
    integer(IK), intent(inout)          :: n_reco
    integer(IK), intent(inout)          :: n_core
    integer(IK), intent(inout)          :: par(n_target)
    integer(IK), intent(inout)          :: sptr(n_heap)
    type(heap_container), intent(inout) :: core(n_edge)
    integer(IK)                         :: i, n_pack
    logical                             :: same
    n_pack = n_reco
    n_reco = 0
    i = 0
    do while (i < n_pack)
      if (hcon(sptr(i + 1))%s /= IS_FINISHED) exit
      i = i + 1
      associate (c => hcon(sptr(i)))
        if (c%lb >= lambda) return
        call unite(n_target, c%p, par, same)
        if (.not. same) then
          n_core = n_core + 1
          core(n_core) = c
        end if
      end associate
      if (n_core == n_edge) return
    end do
    do while (i < n_pack)
      i = i + 1
      associate (c => hcon(sptr(i)))
        if (c%lb >= lambda) return
        if (is_same(n_target, c%p, par)) cycle
        n_reco = n_reco + 1
        sptr(n_reco) = sptr(i)
      end associate
    end do
  end subroutine kruscal
!
  pure subroutine update_heap(n_reco, n_heap, hcon, sptr, hptr)
    integer(IK), intent(in)          :: n_reco
    integer(IK), intent(in)          :: n_heap
    type(heap_container), intent(in) :: hcon(n_heap)
    integer(IK), intent(in)          :: sptr(n_reco)
    integer(IK), intent(inout)       :: hptr(n_reco)
    integer(IK)                      :: i, j
    j = n_reco
    do i = 1, n_reco
      hptr(j) = sptr(i)
      call down_heap(n_reco, j, hcon, hptr)
      j = j - 1
    end do
  end subroutine update_heap
!
  pure elemental subroutine swap(a, b)
    integer(IK), intent(inout) :: a, b
    integer(IK)                :: s
    s = b
    b = a
    a = s
  end subroutine swap
!
! union find
!
  pure subroutine init_arange(n, a)
    integer(IK), intent(in)    :: n
    integer(IK), intent(inout) :: a(n)
    integer(IK)                :: i
    do concurrent(i=1:n)
      a(i) = i
    end do
  end subroutine init_arange
!
  pure subroutine findroot(n, v, par, root, rank)
    integer(IK), intent(in)    :: n, v, par(n)
    integer(IK), intent(inout) :: root, rank
    root = v
    rank = 0
    do
      if (par(root) == root) exit
      root = par(root)
      rank = rank + 1
    end do
  end subroutine findroot
!
  pure subroutine reduction(n, v, root, par)
    integer(IK), intent(in)    :: n, v, root
    integer(IK), intent(inout) :: par(n)
    integer(IK)                :: p
    p = v
    do
      if (par(p) == root) exit
      p = par(p)
      par(p) = root
    end do
  end subroutine reduction
!
  pure function is_same(n, p, par)
    integer(IK), intent(in) :: n, p
    integer(IK), intent(in) :: par(n)
    logical                 :: is_same
    integer(IK)             :: v1, v2, r1, l1, r2, l2
    call cantor_pair_inverse(p, v1, v2)
    call findroot(n, v1, par, r1, l1)
    call findroot(n, v2, par, r2, l2)
    is_same = r1 == r2
  end function is_same
!
  pure subroutine unite(n, p, par, is_same)
    integer(IK), intent(in)    :: n, p
    integer(IK), intent(inout) :: par(n)
    logical, intent(inout)     :: is_same
    integer(IK)                :: v1, v2, r1, l1, r2, l2
    call cantor_pair_inverse(p, v1, v2)
    call findroot(n, v1, par, r1, l1)
    call findroot(n, v2, par, r2, l2)
    is_same = r1 == r2
    if (is_same) then
      call reduction(n, v1, r1, par)
      call reduction(n, v2, r2, par)
    else
      if (l1 > l2) then
        par(r2) = r1
        call reduction(n, v1, r1, par)
        call reduction(n, v2, r1, par)
      else
        par(r1) = r2
        call reduction(n, v2, r2, par)
        call reduction(n, v1, r2, par)
      end if
    end if
  end subroutine unite
!
! pure elemental subroutine cantor_pair(i, j, k)
!   integer(IK), intent(in)    :: i, j
!   integer(IK), intent(inout) :: k
!   if (i > j) then
!     k = (i - 2) * (i - 1) / 2 + j
!   else
!     k = (j - 2) * (j - 1) / 2 + i
!   end if
! end subroutine cantor_pair
!
  pure subroutine make_heap(n, e, h)
    integer(IK), intent(in)          :: n
    type(heap_container), intent(in) :: e(*)
    integer(IK), intent(inout)       :: h(n)
    integer(IK)                      :: i
    do i = n, 1, -1
      call down_heap(n, i, e, h)
    end do
  end subroutine make_heap
!
  pure subroutine down_heap(n, p, e, h)
    integer(IK), intent(in)          :: n, p
    type(heap_container), intent(in) :: e(*)
    integer(IK), intent(inout)       :: h(*)
    integer(IK)                      :: r, s, t
    r = p
    do while (r + r <= n)
      s = r + r
      t = s + 1
      if (t > n) then
        if (e(h(r))%lb > e(h(s))%lb) return
      else
        if (e(h(s))%lb > e(h(t))%lb) then
          if (e(h(r))%lb > e(h(s))%lb) return
        else
          if (e(h(r))%lb > e(h(t))%lb) return
          s = t
        end if
      end if
      call swap(h(r), h(s))
      r = s
    end do
  end subroutine down_heap
!
  pure elemental subroutine cantor_pair_inverse(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    j = (INT(SQRT(real(8 * k - 7))) - 1) / 2 + 2 ! cantor_pair_inverse
    i = k - (j - 1) * (j - 2) / 2
  end subroutine cantor_pair_inverse
!
end module mod_mobbrmsd_mst

