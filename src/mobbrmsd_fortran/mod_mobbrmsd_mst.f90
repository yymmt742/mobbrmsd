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
  implicit none
  private
  public :: mobbrmsd_min_span_tree
!
  type edge_data
    sequence
    integer(IK) :: p, w
    real(RK)    :: lb, ub
    type(mobbrmsd_state) :: state
  end type edge_data
!
! type hash_table
!   sequence
!   integer(IK)          :: k
!   type(mobbrmsd_state) :: state
! end type hash_table
!
! type ij
!   sequence
!   integer(IK) :: i, j
! end type ij
!
contains
!| minimum spanning tree construction.
  subroutine mobbrmsd_min_span_tree( &
 &             n_target, &
 &             header, &
 &             X, &
 &             W, &
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
    real(RK), intent(inout)             :: W(*)
    !! work memory, must be larger than header%memsize() * n_target * (n_target-1) / 2
    logical, intent(in), optional       :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional       :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
    integer(IK), intent(out), optional  :: edges(2, n_target - 1)
    !! minimum spanning tree edges
    real(RK), intent(out), optional     :: weights(n_target - 1)
    !! minimum spanning tree weights
    type(edge_data)                     :: core(n_target - 1)
    integer(IK)                         :: memsize, ldxsize
    integer(IK)                         :: n_edges, n_chunk, n_hash
    integer(IK)                         :: i, j, k
    if (.not. PRESENT(edges) .and. .not. PRESENT(weights)) return
    n_chunk = MIN(n_target * 30, n_target * (n_target - 1) / 2)
    n_hash = n_chunk + 2
    n_edges = n_target - 1
    memsize = mobbrmsd_memsize(header)
    ldxsize = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
    block
      type(edge_data)  :: edge(n_chunk)
      real(RK)         :: lambda
      real(RK)         :: w(n_chunk * memsize)
      integer(IK)      :: heap(n_chunk) ! for heap sort
      integer(IK)      :: par(n_target) ! for union find
      integer(IK)      :: n_core, n_edge, n_pack
      lambda = 1.0E-2_RK
      call init_par(n_target, par)
      call init_edge(edge)
      n_core = 0
      n_edge = 0
      do
        call construct_chunk( &
           &  memsize &
           &, ldxsize &
           &, n_target &
           &, n_edge &
           &, n_core &
           &, n_chunk &
           &, header &
           &, lambda &
           &, par &
           &, X &
           &, W &
           &, n_pack &
           &, heap &
           &, edge &
           &, remove_com &
           &, sort_by_g &
           & )
        block
          integer(IK) :: g
          real(RK) :: lambda1
          g = pick_heap(n_pack, n_edges - n_core, edge, heap)
          lambda1 = MAX(lambda, edge(g)%ub)
          call kruscal( &
              &  n_target &
              &, n_chunk &
              &, n_pack &
              &, lambda &
              &, n_edge &
              &, n_core &
              &, par &
              &, heap &
              &, edge &
              &, core &
              & )
          if (n_core == n_edges) exit
          lambda = lambda1
        end block
      end do
    end block
!
    if (PRESENT(edges)) then
      do concurrent(k=1:n_edges)
        call cantor_pair_inverse(core(k)%p, i, j)
        edges(1, k) = i
        edges(2, k) = j
      end do
    end if
    if (PRESENT(weights)) then
      do concurrent(k=1:n_edges)
        weights(k) = core(k)%lb
      end do
    end if
  end subroutine mobbrmsd_min_span_tree
!
  subroutine construct_chunk( &
            &  memsize, &
            &  ldxsize, &
            &  n_target, &
            &  n_edge, &
            &  n_core, &
            &  n_chunk, &
            &  header, &
            &  lambda, &
            &  par, &
            &  X, &
            &  W, &
            &  n_pack, &
            &  heap, &
            &  edge, &
            &  remove_com, &
            &  sort_by_g &
            & )
    integer(IK), intent(in)        :: memsize
    integer(IK), intent(in)        :: ldxsize
    integer(IK), intent(in)        :: n_target
    integer(IK), intent(in)        :: n_edge
    integer(IK), intent(in)        :: n_core
    integer(IK), intent(in)        :: n_chunk
    integer(IK), intent(in)        :: par(n_target)
    type(mobbrmsd), intent(in)     :: header
    real(RK), intent(in)           :: lambda
    real(RK), intent(in)           :: X(*)
    real(RK), intent(inout)        :: W(*)
    integer(IK), intent(inout)     :: n_pack
    integer(IK), intent(inout)     :: heap(n_chunk)
    type(edge_data), intent(inout) :: edge(n_chunk)
    logical, intent(in), optional  :: remove_com, sort_by_g
    integer(IK)                    :: edges(n_edge + 1)
    logical                        :: check1, check2, check3
    integer(IK)                    :: k, l, addr, ld, lc
!
    ld = n_target * (n_target - 1) / 2
    lc = MIN(n_chunk, ld - n_core)
!
    !$omp parallel do
    do k = 1, n_edge
      call update( &
          &  ldxsize, &
          &  n_chunk, &
          &  header, &
          &  lambda, &
          &  X, &
          &  W, &
          &  edge(heap(k)), &
          &  remove_com, &
          &  sort_by_g &
          &)
    end do
    !$omp end parallel do
    if (n_edge == n_chunk) call make_heap(n_edge, edge, heap)
!
    do concurrent(k=1:n_edge)
      edges(k) = edge(heap(k))%p
    end do
    call sort(n_edge, edges)
    edges(n_edge + 1) = ld + 1
!
    n_pack = n_edge
    l = 0
    k = 1
    !$omp parallel shared(edge,heap,n_pack,k,l), private(addr,check1,check2)
    do
      !$omp critical
      l = l + 1
      check1 = n_pack >= n_chunk
      check2 = l > ld
      if (check2) then
        check3 = .false.
      else
        check3 = l == edges(k)
        if (check3) then
          k = k + 1
        else
          check3 = is_same(n_target, l, par)
        end if
      end if
      if (.not. (check1 .or. check2 .or. check3)) then
        n_pack = n_pack + 1
        addr = find_address(n_chunk, l, edge)
        heap(n_pack) = addr
        edge(addr)%p = -l
        edge(addr)%w = (addr - 1) * memsize + 1
      else
        addr = 0
      end if
      !$omp end critical

      if (check1) exit
      if (check2) exit
      if (check3) cycle

      call update( &
          &  ldxsize, &
          &  n_chunk, &
          &  header, &
          &  lambda, &
          &  X, &
          &  W, &
          &  edge(addr), &
          &  remove_com, &
          &  sort_by_g &
          &)
    end do
    !$omp end parallel

    if (n_pack == n_chunk) call make_heap(n_pack, edge, heap)

    !$omp parallel shared(edge,heap,k,l), private(check1,check2)
    do
      block
        real(RK) :: V(memsize)
        type(edge_data) :: t
        !$omp critical
        l = l + 1
        check1 = l > ld
        if (check1) then
          check2 = .false.
        else
          check2 = l == edges(k)
          if (check2) then
            k = k + 1
          else
            check2 = is_same(n_target, l, par)
          end if
        end if
        if (.not. (check1 .or. check2)) then
          call init_edge(t)
          t%p = -l
        end if
        !$omp end critical

        if (check1) exit
        if (check2) cycle

        call update( &
            &  ldxsize, &
            &  n_chunk, &
            &  header, &
            &  lambda, &
            &  X, &
            &  V, &
            &  t, &
            &  remove_com, &
            &  sort_by_g &
            &)
        !$omp critical
        if (t%lb < edge(heap(1))%lb) then
          edge(heap(1)) = t
          edge(heap(1))%w = (heap(1) - 1) * memsize + 1
          call memcopy(memsize, V, W(edge(heap(1))%w))
          call down_heap(n_chunk, 1, edge, heap)
        end if
        !$omp end critical
      end block
    end do
    !$omp end parallel

    if (n_pack < n_chunk) call make_heap(n_pack, edge, heap)

!   block
!     real(RK) :: lambda_local
!     integer(IK) :: g, n_rem
!     n_rem = n_target - 1 - n_core
!     g = pick_heap(n_pack, n_rem, edge, heap)
!     lambda_local = edge(g)%ub
!     do while (lambda_local < lambda .and. .not. mobbrmsd_state_is_finished(edge(g)%state))
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
  end subroutine construct_chunk
!
  pure subroutine update( &
 &             ldxsize, &
 &             n_chunk, &
 &             header, &
 &             lambda, &
 &             X, &
 &             W, &
 &             edge, &
 &             remove_com, &
 &             sort_by_g &
 &            )
    integer(IK), intent(in)        :: ldxsize
    integer(IK), intent(in)        :: n_chunk
    type(mobbrmsd), intent(in)     :: header
    real(RK), intent(in)           :: lambda, X(*)
    real(RK), intent(inout)        :: W(*)
    type(edge_data), intent(inout) :: edge
    logical, intent(in), optional  :: remove_com, sort_by_g
    if (edge%p < 0) then
      block
        integer(IK)          :: i, j, xpnt, ypnt
        edge%p = -edge%p
        call cantor_pair_inverse(edge%p, i, j)
        xpnt = (i - 1) * ldxsize + 1
        ypnt = (j - 1) * ldxsize + 1
        call mobbrmsd_run( &
       &       header, &
       &       edge%state,  &
       &       X(xpnt), &
       &       X(ypnt), &
       &       W(edge%w), &
       &       cutoff=lambda, &
       &       sort_by_g=sort_by_g &
       &      )
      end block
    else
      call mobbrmsd_restart( &
     &       header, &
     &       edge%state,  &
     &       W(edge%w), &
     &       cutoff=lambda &
     &      )
    end if
    edge%lb = mobbrmsd_state_lowerbound_as_rmsd(edge%state)
    edge%ub = mobbrmsd_state_upperbound_as_rmsd(edge%state)
  end subroutine update
!
  pure function pick_heap(n, p, e, h) result(res)
    integer(IK), intent(in)     :: n, p
    type(edge_data), intent(in) :: e(n)
    integer(IK), intent(in)     :: h(n)
    integer(IK)                 :: t(n)
    integer(IK)                 :: k, res
    do concurrent(k=1:n)
      t(k) = h(k)
    end do
    do k = n, n - p + 1, -1
      t(1) = h(k)
      call down_heap(k - 1, 1, e, t)
    end do
    res = t(1)
  end function pick_heap
!
  pure subroutine kruscal(n_target, n_chunk, n_pack, lambda, n_edge, n_core, par, heap, edge, core)
    integer(IK), intent(in)        :: n_target, n_chunk, n_pack
    real(RK), intent(in)           :: lambda
    integer(IK), intent(inout)     :: n_edge, n_core
    integer(IK), intent(inout)     :: par(n_target)
    integer(IK), intent(inout)     :: heap(n_chunk)
    type(edge_data), intent(inout) :: edge(n_chunk), core(n_target - 1)
    integer(IK)                    :: t(n_pack)
    integer(IK)                    :: k
    logical                        :: is_same
    do concurrent(k=1:n_pack)
      t(k) = heap(k)
    end do
    do concurrent(k=1:n_chunk)
      heap(k) = 0
    end do
    do k = n_pack, 2, -1
      call swap(t(1), t(k))
      call down_heap(k - 1, 1, edge, t)
    end do
    n_edge = 0
    do k = 1, n_pack
      if (edge(t(k))%ub < lambda) then
        call unite(n_target, edge(t(k))%p, par, is_same)
        if (is_same) then
          call init_edge(edge(t(k)))
        else
          n_core = n_core + 1
          core(n_core) = edge(t(k))
          call init_edge(edge(t(k)))
          if (n_core == n_target - 1) exit
        end if
      else
        n_edge = n_edge + 1
        heap(n_edge) = t(k)
      end if
    end do
  end subroutine kruscal
!
  pure elemental subroutine swap(a, b)
    integer(IK), intent(inout) :: a, b
    integer(IK)                :: s
    s = b
    b = a
    a = s
  end subroutine swap
!
  pure subroutine down_heap(n, p, e, h)
    integer(IK), intent(in)     :: n, p
    type(edge_data), intent(in) :: e(*)
    integer(IK), intent(inout)  :: h(*)
    integer(IK)                 :: r, s, t
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
! pure recursive subroutine insert(p, q, a)
!   integer, intent(in)       :: p, q
!   type(edge_data), intent(inout) :: a(*)
!   integer                   :: r
!   if (a(p)%lb < a(q)%lb) then
!     if (a(q)%l < 1) then
!       a(q)%l = p
!       return
!     end if
!     r = a(q)%l
!   else
!     if (a(q)%r < 1) then
!       a(q)%r = p
!       return
!     end if
!     r = a(q)%r
!   end if
!   call insert(p, r, a)
! end subroutine insert
!
! pure recursive subroutine ordering(p, q, e)
!   integer(IK), intent(in)    :: p
!   integer(IK), intent(inout) :: q
!   type(edge_data), intent(inout)  :: e(*)
!   if (e(p)%r > 0) call ordering(e(p)%r, q, e)
!   e(p)%s = q
!   q = p
!   if (e(p)%l > 0) call ordering(e(p)%l, q, e)
! end subroutine ordering
!
! union find
!
  pure subroutine init_par(n, par)
    integer(IK), intent(in)    :: n
    integer(IK), intent(inout) :: par(n)
    integer(IK)                :: i
    do concurrent(i=1:n)
      par(i) = i
    end do
  end subroutine init_par
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
  pure function find_address(n, p, e) result(res)
    integer(IK), intent(in)     :: n, p
    type(edge_data), intent(in) :: e(n)
    integer(IK)                 :: c, i, res
    res = 0
    i = hash(n, p)
    c = i
    do
      if (e(c)%p == p) then
        res = c
        return
      elseif (e(c)%p == 0) then
        if (res == 0) res = c
      end if
      c = MODULO(c, n) + 1
      if (c == i) return
    end do
  end function find_address

! pure elemental subroutine hash_table_init(h)
!   type(hash_table), intent(inout) :: h
!   h%k = 0
! end subroutine hash_table_init

! pure subroutine hash_add(n, k, h)
!   integer(IK), intent(in)         :: n, k
!   real(RK), intent(in)            :: v
!   type(hash_table), intent(inout) :: h(n)
!   integer(IK)                     :: c
!   c = hash(n, k)
!   do
!     if (h(c)%k == k) return
!     if (h(c)%k <= 0) then
!       h(c)%k = k
!       return
!     else
!       c = MODULO(c, n) + 1
!     end if
!   end do
! end subroutine hash_add

! pure function hash_get(n, k, h) result(c)
!   integer(IK), intent(in)      :: n, k
!   type(hash_table), intent(in) :: h(n)
!   integer(IK)                  :: c, i
!   i = hash(n, k)
!   c = i
!   do
!     if (h(c)%k == k) then
!       return
!     elseif (h(c)%k == 0) then
!       c = 0
!       return
!     end if
!     c = MODULO(c, n) + 1
!     if (c == i) then
!       c = 0
!       return
!     end if
!   end do
! end function hash_get

! pure subroutine hash_replace(n, k, l, h)
!   integer(IK), intent(in)         :: n, k, l
!   type(hash_table), intent(inout) :: h(n)
!   integer(IK)                     :: c, i
!   i = hash(n, k)
!   c = i
!   do
!     if (h(c)%k == k) then
!       h(c)%k = l
!       return
!     elseif (h(c)%k == 0) then
!       h(c)%k = l
!       return
!     end if
!     c = MODULO(c, n) + 1
!     if (c == i) return
!   end do
! end subroutine hash_replace

! pure subroutine hash_del(n, k, h)
!   integer(IK), intent(in)         :: n, k
!   type(hash_table), intent(inout) :: h(n)
!   integer(IK)                     :: c, i
!   i = hash(n, k)
!   c = i
!   do
!     if (h(c)%k == k) then
!       h(c)%k = -k
!       return
!     elseif (h(c)%k == 0) then
!       return
!     end if
!     c = MODULO(c, n) + 1
!     if (c == i) return
!   end do
! end subroutine hash_del
!
! cantor_pair
!
  pure subroutine make_heap(n, e, h)
    integer(IK), intent(in)     :: n
    type(edge_data), intent(in) :: e(n)
    integer(IK), intent(inout)  :: h(n)
    integer(IK)                 :: i
    do concurrent(i=1:n)
      h(i) = i
    end do
    do i = n, 1, -1
      call down_heap(n, i, e, h)
    end do
  end subroutine make_heap
!
  pure elemental subroutine cantor_pair_inverse(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    j = (INT(SQRT(real(8 * k - 7))) - 1) / 2 + 2 ! cantor_pair_inverse
    i = k - (j - 1) * (j - 2) / 2
  end subroutine cantor_pair_inverse
!
  pure elemental subroutine init_edge(e)
    type(edge_data), intent(inout) :: e
    e%p = 0
    e%w = 1
    e%lb = HUGE(0.0_RK)
    e%ub = HUGE(0.0_RK)
  end subroutine init_edge
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

