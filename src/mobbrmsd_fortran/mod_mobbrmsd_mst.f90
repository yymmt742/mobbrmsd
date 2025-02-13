!| Configure a minimum global tree with mobbrmsd. <br>
!  MST construction is based on the Prims algorithm,
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
  type heap_container
    sequence
    integer(IK) :: p, s
    real(RK)    :: lb, ub
  end type heap_container
!
contains
!| minimum spanning tree construction.
  subroutine mobbrmsd_min_span_tree( &
 &             n_target, &
 &             header, &
 &             X, &
 &             n_chunk, &
 &             remove_com, &
 &             sort_by_g, &
 &             edges, &
 &             weights &
 &            )
    integer(IK), intent(in)             :: n_target
    !! number of coordinates
    type(mobbrmsd), intent(in)          :: header
    !! mobbrmsd_header
    real(RK), intent(in)                :: X(*)
    !! coordinate sequence
    integer(IK), intent(in), optional   :: n_chunk
    !! Maximum memory storage area size
    logical, intent(in), optional       :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional       :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
    integer(IK), intent(out), optional  :: edges(2, n_target - 1)
    !! minimum spanning tree edges
    real(RK), intent(out), optional     :: weights(n_target - 1)
    !! minimum spanning tree weights
    integer(IK)                         :: par(n_target) ! for union find
    integer(IK)                         :: n_tri, ldw, ldx
    integer(IK)                         :: n_core, n_reco, n_thread, n_edge, n_heap
    integer(IK)                         :: i, j, k
!
    if (.not. PRESENT(edges) .and. .not. PRESENT(weights)) return
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
    if (PRESENT(n_chunk)) then
      if (n_chunk > 0) n_heap = MIN(MAX(n_target * 2, n_chunk), n_heap)
    end if
    n_heap = MAX(n_heap, n_thread)
!
    block
      type(mobbrmsd_state), allocatable :: state(:)
      integer(IK)          :: hptr(n_heap) ! heap pointer
      integer(IK)          :: sorted(n_edge) ! top + core
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
                               &, maxeval=1 &
                               &, remove_com=remove_com &
                               &, sort_by_g=sort_by_g &
                               & )
!
      n_core = 0
      n_reco = 0
      lambda = 0.5_RK
      call init_arange(n_target, par)
      call init_arange(n_edge, sorted)
      call init_arange(n_heap, hptr)
!
      do while (n_core < n_edge)
!
        call reflesh_heap(n_tri &
                       &, n_target &
                       &, n_edge &
                       &, n_heap &
                       &, n_core &
                       &, par &
                       &, state &
                       &, n_reco &
                       &, hptr &
                       &, hcon &
                       &)
!
        do while (n_core < n_edge)
          call sort_heap(n_tri &
                      &, n_edge &
                      &, n_heap &
                      &, n_reco &
                      &, n_core &
                      &, state &
                      &, hptr &
                      &, hcon &
                      &, sorted &
                      & )
          print'(*(f8.4))', hcon(sorted)%ub
          print'(*(f8.4))', hcon(sorted)%lb
          stop
!
!         call kruscal(n_target &
!                   &, n_heap &
!                   &, n_pack &
!                   &, lambda &
!                   &, n_reco &
!                   &, n_core &
!                   &, par &
!                   &, heap &
!                   &, edge &
!                   &, core &
!                   &)
          call update_bounds(ldx &
                          &, ldw &
                          &, n_edge &
                          &, n_heap &
                          &, header &
                          &, lambda &
                          &, X &
                          &, wrkmem &
                          &, sorted &
                          &, hcon &
                          &, state &
                          &, remove_com &
                          &, sort_by_g &
                          &)
!
        end do
      end do
!
!   if (PRESENT(edges)) then
!     do concurrent(k=1:n_edge)
!       call cantor_pair_inverse(core(k)%p, i, j)
!       edges(1, k) = i
!       edges(2, k) = j
!     end do
!   end if
!   if (PRESENT(weights)) then
!     do concurrent(k=1:n_edge)
!       weights(k) = core(k)%lb
!     end do
!   end if
    end block
  end subroutine mobbrmsd_min_span_tree
!
  pure subroutine reflesh_heap(n_tri &
                            &, n_target &
                            &, n_edge &
                            &, n_heap &
                            &, n_core &
                            &, par &
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
    integer(IK), intent(in)             :: par(n_edge)
    type(mobbrmsd_state), intent(in)    :: state(n_tri)
    integer(IK), intent(inout)          :: n_reco
    integer(IK), intent(inout)          :: hptr(n_heap)
    type(heap_container), intent(inout) :: hcon(n_heap)
    real(RK)                            :: lb
    integer(IK)                         :: record(n_reco + 1)
    integer(IK)                         :: i, j, n_rem
!
    n_rem = MIN(n_heap, n_tri - n_core)
    do concurrent(i=1:n_reco)
      record(i) = hcon(hptr(i))%p
    end do
    call sort(n_reco, record)
    record(n_reco + 1) = n_tri + 1
!
    i = 0
    j = 1
    do while (i < n_tri .and. n_reco < n_heap)
      i = i + 1
      if (i == record(j)) then
        j = j + 1
        cycle
      end if
      if (is_same(n_target, i, par)) cycle
      n_reco = n_reco + 1
      hcon(hptr(n_reco))%p = i
      hcon(hptr(n_reco))%s = MERGE(-1, 0, mobbrmsd_state_is_finished(state(i)))
      hcon(hptr(n_reco))%lb = mobbrmsd_state_lowerbound_as_rmsd(state(i))
      hcon(hptr(n_reco))%ub = mobbrmsd_state_upperbound_as_rmsd(state(i))
    end do
    call make_heap(n_reco, hcon, hptr)
    do while (i < n_tri)
      if (is_same(n_target, i, par)) cycle
      lb = mobbrmsd_state_lowerbound_as_rmsd(state(i))
      if (lb < hcon(hptr(1))%lb) then
        hcon(hptr(1))%p = i
        hcon(hptr(1))%s = MERGE(-1, 0, mobbrmsd_state_is_finished(state(i)))
        hcon(hptr(1))%lb = lb
        hcon(hptr(1))%ub = mobbrmsd_state_lowerbound_as_rmsd(state(i))
        call down_heap(n_heap, 1, hcon, hptr)
      end if
    end do
  end subroutine reflesh_heap
!
  subroutine sort_heap(n_tri &
                    &, n_edge &
                    &, n_heap &
                    &, n_reco &
                    &, n_core &
                    &, state &
                    &, hptr &
                    &, hcon &
                    &, sorted &
                    & )
    integer(IK), intent(in)          :: n_tri
    integer(IK), intent(in)          :: n_edge
    integer(IK), intent(in)          :: n_heap
    integer(IK), intent(in)          :: n_reco
    integer(IK), intent(in)          :: n_core
    type(mobbrmsd_state), intent(in) :: state(n_tri)
    integer(IK), intent(in)          :: hptr(n_reco)
    type(heap_container), intent(in) :: hcon(n_heap)
    integer(IK), intent(inout)       :: sorted(n_edge)
    integer(IK)                      :: t(n_heap)
    integer(IK)                      :: k
    do concurrent(k=1:n_reco)
      t(k) = hptr(k)
    end do
    do k = n_reco, 2, -1
      call swap(t(1), t(k))
      call down_heap(k - 1, 1, hcon, t)
    end do
    do concurrent(k=1:n_edge - n_core)
      sorted(k) = t(k)
    end do
  end subroutine sort_heap
!
  subroutine update_bounds(ldx &
                        &, ldw &
                        &, n_edge &
                        &, n_heap &
                        &, header &
                        &, lambda &
                        &, X &
                        &, wrkmem &
                        &, hptr &
                        &, hcon &
                        &, state &
                        &, remove_com &
                        &, sort_by_g &
                        & )
    integer(IK), intent(in)             :: ldx
    integer(IK), intent(in)             :: ldw
    integer(IK), intent(in)             :: n_edge
    integer(IK), intent(in)             :: n_heap
    type(mobbrmsd), intent(in)          :: header
    real(RK), intent(in)                :: lambda
    real(RK), intent(in)                :: X(*)
    real(RK), intent(inout)             :: wrkmem(ldw, n_heap)
    integer(IK), intent(in)             :: hptr(n_heap)
    type(heap_container), intent(inout) :: hcon(n_heap)
    type(mobbrmsd_state), intent(inout) :: state(*)
    logical, intent(in), optional       :: remove_com, sort_by_g
    integer(IK)                         :: k
    !$omp parallel do
    do k = 1, n_edge
      associate (p => hcon(hptr(k))%p&
              &, s => hcon(hptr(k))%s&
              &, h => hptr(k)&
              & )
      if (s == 0) then
        block
          integer(IK) :: i, j
          call cantor_pair_inverse(p, i, j)
          call mobbrmsd_run(header &
                         &, state(p) &
                         &, X((i - 1) * ldx + 1) &
                         &, X((j - 1) * ldx + 1) &
                         &, wrkmem(1, h) &
                         &, cutoff=lambda &
                         &, remove_com=remove_com &
                         &, sort_by_g=sort_by_g &
                         & )
        end block
      elseif (s == 1) then
        call mobbrmsd_restart(header &
                           &, state(p) &
                           &, wrkmem(1, h) &
                           &, cutoff=lambda &
                           & )
      end if
      hcon(h)%lb = mobbrmsd_state_lowerbound_as_rmsd(state(p))
      hcon(h)%ub = mobbrmsd_state_lowerbound_as_rmsd(state(p))
      s = MERGE(-1, 1, mobbrmsd_state_is_finished(state(p)))
      end associate
    end do
    !$omp end parallel do
  end subroutine update_bounds
!
! subroutine construct_chunk(ld &
!                         &, memsize &
!                         &, ldxsize &
!                         &, n_target &
!                         &, n_rec &
!                         &, n_core &
!                         &, n_chunk &
!                         &, header &
!                         &, lambda &
!                         &, par &
!                         &, first &
!                         &, X &
!                         &, W &
!                         &, n_pack &
!                         &, heap &
!                         &, edge &
!                         &, state &
!                         &, remove_com &
!                         &, sort_by_g &
!                         & )
!   integer(IK), intent(in)             :: ld
!   integer(IK), intent(in)             :: memsize
!   integer(IK), intent(in)             :: ldxsize
!   integer(IK), intent(in)             :: n_target
!   integer(IK), intent(in)             :: n_rec
!   integer(IK), intent(in)             :: n_core
!   integer(IK), intent(in)             :: n_chunk
!   integer(IK), intent(in)             :: par(n_target)
!   logical, intent(in)                 :: first
!   type(mobbrmsd), intent(in)          :: header
!   real(RK), intent(in)                :: lambda
!   real(RK), intent(in)                :: X(*)
!   real(RK), intent(inout)             :: W(*)
!   integer(IK), intent(inout)          :: n_pack
!   integer(IK), intent(inout)          :: heap(n_chunk)
!   type(edge_data), intent(inout)      :: edge(n_chunk)
!   type(mobbrmsd_state), intent(inout) :: state(*)
!   logical, intent(in), optional       :: remove_com, sort_by_g
!   integer(IK)                         :: edges(n_rec + 1)
!   logical                             :: check1, check2, check3
!   integer(IK)                         :: k, l, addr, n_rem
!
!   n_rem = MIN(n_chunk, ld - n_core)
!   do concurrent(k=1:n_rec)
!     edges(k) = edge(heap(k))%p
!   end do
!   call sort(n_rec, edges)
!   edges(n_rec + 1) = ld + 1
!
!   !$omp parallel do
!   do k = 1, n_rec
!     call update(ldxsize &
!              &, n_chunk &
!              &, header &
!              &, lambda &
!              &, first &
!              &, X &
!              &, W &
!              &, state &
!              &, edge(heap(k)) &
!              &, remove_com &
!              &, sort_by_g &
!              &)
!   end do
!   !$omp end parallel do
!   n_pack = n_rec
!   l = 0
!   k = 1
!
!   if (n_rec < n_rem) then
!
!     !$omp parallel shared(edge,heap,n_pack,k,l), private(addr,check1,check2,check3)
!     do
!       !$omp critical
!       if (l <= ld) l = l + 1
!       check1 = n_pack >= n_rem
!       check2 = l > ld
!       if (check2) then
!         check3 = .false.
!       else
!         check3 = l == edges(k)
!         if (check3) then
!           k = k + 1
!         else
!           check3 = is_same(n_target, l, par)
!         end if
!       end if
!       if (check1 .or. check2 .or. check3) then
!         addr = 0
!       else
!         n_pack = n_pack + 1
!         addr = find_address(n_chunk, l, edge)
!         heap(n_pack) = addr
!         edge(addr)%p = -l
!         !edge(addr)%w = (addr - 1) * memsize + 1
!       end if
!       !$omp end critical

!       !print *, l, check1, check2, check3
!       if (check1) exit
!       if (check2) exit
!       if (check3) cycle

!       call update(ldxsize &
!                &, n_chunk &
!                &, header &
!                &, lambda &
!                &, first &
!                &, X &
!                &, W &
!                &, state &
!                &, edge(addr) &
!                &, remove_com &
!                &, sort_by_g &
!                &)
!     end do
!     !$omp end parallel
!   end if

!   call make_heap(n_pack, edge, heap)

!   !$omp parallel shared(edge,heap,k,l), private(check1,check2)
!   do
!     block
!       real(RK) :: V(memsize)
!       type(edge_data) :: t
!       !$omp critical
!       if (l <= ld) l = l + 1
!       check1 = l > ld
!       if (check1) then
!         check2 = .false.
!       else
!         check2 = l == edges(k)
!         if (check2) then
!           k = k + 1
!         else
!           check2 = is_same(n_target, l, par)
!         end if
!       end if
!       if (.not. (check1 .or. check2)) then
!         call init_edge(t)
!         t%p = -l
!       end if
!       !$omp end critical

!       if (check1) exit
!       if (check2) cycle

!       call update(ldxsize &
!                &, n_chunk &
!                &, header &
!                &, lambda &
!                &, first &
!                &, X &
!                &, V &
!                &, state &
!                &, t &
!                &, remove_com &
!                &, sort_by_g &
!                &)
!       !$omp critical
!       if (t%lb < edge(heap(1))%lb) then
!         t%w = edge(heap(1))%w
!         edge(heap(1)) = t
!         !edge(heap(1))%w = (heap(1) - 1) * memsize + 1
!         call memcopy(memsize, V, W(edge(heap(1))%w))
!         call down_heap(n_chunk, 1, edge, heap)
!       end if
!       !$omp end critical
!     end block
!   end do
!   !$omp end parallel

!   block
!     real(RK) :: lambda_local
!     integer(IK) :: g, n_rem
!     n_rem = n_target - 1 - n_core
!     g = pick_heap(n_pack, n_rem, edge, heap)
!     lambda_local = edge(g)%ub
!     do while (lambda_local < lambda .and. mst_is_unfinished(n_pack, n_rem, heap, edge))
!       do k = 1, n_rem
!         if (mobbrmsd_state_is_finished(edge(heap(k))%state)) cycle
!         call update( &
!             &  ldxsize, &
!             &  n_chunk, &
!             &  header, &
!             &  lambda_local, &
!             &  X, &
!             &  W, &
!             &  edge(heap(k)), &
!             &  remove_com, &
!             &  sort_by_g &
!             & )
!       end do
!       call make_heap(n_pack, edge, heap)
!       g = pick_heap(n_pack, n_rem, edge, heap)
!       lambda_local = MAX(lambda_local, edge(g)%ub)
!     end do
!   end block
! end subroutine construct_chunk
!
! pure function mst_is_unfinished(n_heap, n_rem, h, e) result(res)
!   integer(IK), intent(in)     :: n_heap, n_rem
!   integer(IK), intent(in)     :: h(n_heap)
!   type(edge_data), intent(in) :: e(*)
!   logical                     :: res
!   integer(IK)                 :: t(n_heap)
!   integer(IK)                 :: k
!   do concurrent(k=1:n_heap)
!     t(k) = h(k)
!   end do
!   do k = n_heap, 2, -1
!     call swap(t(1), t(k))
!     call down_heap(k - 1, 1, e, t)
!   end do
!   res = .not. ALL(mobbrmsd_state_is_finished(e(t(:n_rem))%state))
! end function mst_is_unfinished
!
! subroutine update(ldxsize &
!                &, n_chunk &
!                &, header &
!                &, lambda &
!                &, first &
!                &, X &
!                &, W &
!                &, state &
!                &, edge &
!                &, remove_com &
!                &, sort_by_g &
!                & )
!   integer(IK), intent(in)             :: ldxsize
!   integer(IK), intent(in)             :: n_chunk
!   type(mobbrmsd), intent(in)          :: header
!   logical, intent(in)                 :: first
!   real(RK), intent(in)                :: lambda, X(*)
!   real(RK), intent(inout)             :: W(*)
!   type(mobbrmsd_state), intent(inout) :: state(*)
!   type(heap_container), intent(inout) :: edge
!   logical, intent(in), optional       :: remove_com, sort_by_g
!   if (first) then
!     block
!       integer(IK)          :: i, j, xpnt, ypnt
!       call cantor_pair_inverse(edge%p, i, j)
!       xpnt = (i - 1) * ldxsize + 1
!       ypnt = (j - 1) * ldxsize + 1
!       call mobbrmsd_run( &
!      &       header, &
!      &       state(edge%p), &
!      &       X(xpnt), &
!      &       X(ypnt), &
!      &       W(edge%w), &
!      &       cutoff=lambda, &
!      &       sort_by_g=sort_by_g &
!      &      )
!     end block
!   else
!     call mobbrmsd_restart( &
!    &       header, &
!    &       state(edge%p), &
!    &       W(edge%w), &
!    &       cutoff=lambda &
!    &      )
!   end if
!   !print *, edge%p, edge%w, edge%lb, edge%ub,&
!   !  & mobbrmsd_state_lowerbound_as_rmsd(state(edge%p)), mobbrmsd_state_upperbound_as_rmsd(state(edge%p))
!   edge%lb = mobbrmsd_state_lowerbound_as_rmsd(state(edge%p))
!   edge%ub = mobbrmsd_state_upperbound_as_rmsd(state(edge%p))
!   edge%not_finished = .not. mobbrmsd_state_is_finished(state(edge%p))
! end subroutine update
!
  pure function pick_heap(n, p, e, h) result(res)
    integer(IK), intent(in)          :: n, p
    type(heap_container), intent(in) :: e(*)
    integer(IK), intent(in)          :: h(n)
    integer(IK)                      :: t(n)
    integer(IK)                      :: k, res
    do concurrent(k=1:n)
      t(k) = h(k)
    end do
    do k = n, n - p + 1, -1
      t(1) = t(k)
      call down_heap(k - 1, 1, e, t)
    end do
    res = t(1)
  end function pick_heap
!
! subroutine kruscal(&
!                &  n_target &
!                &, n_chunk &
!                &, n_pack &
!                &, lambda &
!                &, n_rec &
!                &, n_core &
!                &, par &
!                &, heap &
!                &, edge &
!                &, core &
!                &)
!   integer(IK), intent(in)             :: n_target, n_chunk, n_pack
!   real(RK), intent(in)                :: lambda
!   integer(IK), intent(inout)          :: n_rec, n_core
!   integer(IK), intent(inout)          :: par(n_target)
!   integer(IK), intent(inout)          :: heap(n_chunk)
!   type(heap_container), intent(inout) :: edge(n_chunk), core(n_target - 1)
!   integer(IK)                         :: t(n_pack)
!   integer(IK)                         :: k
!   logical                             :: is_same
!   !print'(*(f6.3))', edge(heap(:n_pack))%ub
!   !print'(*(f6.3))', edge(heap(:n_pack))%lb
!   do concurrent(k=1:n_pack)
!     t(k) = heap(k)
!   end do
!   do concurrent(k=1:n_chunk)
!     heap(k) = 0
!   end do
!   do k = n_pack, 2, -1
!     call swap(t(1), t(k))
!     call down_heap(k - 1, 1, edge, t)
!   end do
!   k = 1
!   do while (k <= n_pack .and. n_core < n_target - 1)
!     if (edge(t(k))%ub > lambda .or. edge(t(k))%not_finished) exit
!     !print *, 'krucor', edge(t(k))%ub, lambda, edge(t(k))%not_finished
!     call unite(n_target, edge(t(k))%p, par, is_same)
!     if (.not. is_same) then
!       n_core = n_core + 1
!       core(n_core) = edge(t(k))
!     end if
!     edge(t(k))%p = 0
!     edge(t(k))%lb = HUGE(0.0_RK)
!     edge(t(k))%ub = HUGE(0.0_RK)
!     k = k + 1
!   end do
!   n_rec = 0
!   if (n_core == n_target - 1) return
!   do while (k <= n_pack)
!     n_rec = n_rec + 1
!     !print *, 'krurec', edge(t(k))%ub, lambda, edge(t(k))%not_finished
!     heap(n_rec) = t(k)
!     k = k + 1
!   end do
!   ! print *, n_rec, n_core
!   !print '(*(I6))', edge(heap(:n_rec))%p, 0, core(:n_core)%p
!   !print '(*(I6))', edge(heap(:n_rec))%w, 0, core(:n_core)%w
!   !print'(*(f6.3))', edge(heap(:n_rec))%ub, 0.0_RK, core(:n_core)%ub
!   !print'(*(f6.3))', edge(heap(:n_rec))%lb, 0.0_RK, core(:n_core)%lb
! end subroutine kruscal
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
!!! --- hash table ---
! key value must be k>0.
!
  pure elemental function hash(n, p)
    integer(IK), intent(in)    :: n, p
    integer(IK)                :: hash
    hash = MODULO(p - 1, n) + 1
  end function hash
!
! pure function find_address(n, p, e) result(res)
!   integer(IK), intent(in)     :: n, p
!   type(edge_data), intent(in) :: e(n)
!   integer(IK)                 :: c, i, res
!   res = 0
!   i = hash(n, p)
!   c = i
!   do
!     if (e(c)%p == p) then
!       res = c
!       return
!     elseif (e(c)%p == 0) then
!       if (res == 0) res = c
!     end if
!     c = MODULO(c, n) + 1
!     if (c == i) return
!   end do
! end function find_address
!
  pure elemental subroutine heap_container_init(e)
    type(heap_container), intent(inout) :: e
    e%p = 0
    e%s = 1
    e%lb = HUGE(0.0_RK)
    e%ub = HUGE(0.0_RK)
  end subroutine heap_container_init
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
  pure subroutine sort(n, v)
    integer(IK), intent(in)    :: n
    integer(IK), intent(inout) :: v(n)
    integer(IK)                :: i, j, r, s, t
    do i = n, 1, -1
      r = i
      do while (r + r <= n)
        s = r + r
        t = s + 1
        if (t > n) then
          if (v(r) > v(s)) exit
        else
          if (v(s) > v(t)) then
            if (v(r) > v(s)) exit
          else
            if (v(r) > v(t)) exit
            s = t
          end if
        end if
        j = v(r)
        v(r) = v(s)
        v(s) = j
        r = s
      end do
    end do
    do i = n, 2, -1
      j = v(1)
      v(1) = v(i)
      v(i) = j
      r = 1
      do while (r + r <= i - 1)
        s = r + r
        t = s + 1
        if (t > i - 1) then
          if (v(r) > v(s)) exit
        else
          if (v(s) > v(t)) then
            if (v(r) > v(s)) exit
          else
            if (v(r) > v(t)) exit
            s = t
          end if
        end if
        j = v(r)
        v(r) = v(s)
        v(s) = j
        r = s
      end do
    end do
  end subroutine sort
!
  pure elemental subroutine cantor_pair_inverse(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    j = (INT(SQRT(real(8 * k - 7))) - 1) / 2 + 2 ! cantor_pair_inverse
    i = k - (j - 1) * (j - 2) / 2
  end subroutine cantor_pair_inverse
!
  pure subroutine memcopy(n, from, dest)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: from(n)
    real(RK), intent(inout) :: dest(n)
    integer(IK)             :: i
    do concurrent(i=1:n)
      dest(i) = from(i)
    end do
  end subroutine memcopy
!
end module mod_mobbrmsd_mst

