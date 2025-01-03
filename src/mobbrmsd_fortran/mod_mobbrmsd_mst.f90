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
  type edge
    sequence
    integer(IK) :: p, w
    logical     :: f !, t
    real(RK)    :: lb, ub
  end type edge
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
 &             state, &
 &             X, &
 &             W, &
 &             cutoff, &
 &             ub_cutoff, &
 &             difflim, &
 &             maxeval, &
 &             remove_com, &
 &             sort_by_g, &
 &             difflim_absolute, &
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
    real(RK), intent(in), optional      :: ub_cutoff
    !! The search ends when upperbound is determined to be greater than to ub_cutoff.
    real(RK), intent(in), optional      :: difflim
    !! The search ends when the difference between the lower and upper bounds is less than difflim.
    integer(IK), intent(in), optional   :: maxeval
    !! The search ends when ncount exceeds maxiter.
    logical, intent(in), optional       :: remove_com
    !! if true, remove centroids. default [.true.]
    logical, intent(in), optional       :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate. default [.true.]
    logical, intent(in), optional       :: difflim_absolute
    !! if true, use absolute difflim. default [.false.]
    integer(IK), intent(out), optional  :: edges(2, n_target - 1)
    !! minimum spanning tree edges
    real(RK), intent(out), optional     :: weights(n_target - 1)
    !! minimum spanning tree weights
    type(edge)                          :: core(n_target - 1)
    real(RK)                            :: lambda
    integer(IK)                         :: par(n_target) ! for union find
    integer(IK)                         :: n_core, n_chunk, n_pack, i, j, k
!
    n_core = 0
    call init_par(n_target, par)
    n_chunk = MIN(n_target * 300, n_target * (n_target - 1) / 2)
    lambda = 1.0E-8_RK
    block
      type(edge) :: e(n_chunk), g
      do
        call construct_chunk( &
           &  n_target, &
           &  n_core, &
           &  n_chunk, &
           &  header, &
           &  lambda, &
           &  par, &
           &  X, &
           &  n_pack, &
           &  e, &
           &  remove_com, &
           &  sort_by_g &
           & )
        call update_chunk( &
           &  n_target, &
           &  n_pack, &
           &  n_core, &
           &  n_chunk, &
           &  header, &
           &  lambda, &
           &  par, &
           &  X, &
           &  e, &
           &  remove_com, &
           &  sort_by_g &
           & )
        call kruscal_core(n_target, n_pack, n_core, par, e, core)
        if (n_core == n_target - 1) exit
        g = pick_heap(n_pack, n_target - n_core - 1, e)
        lambda = MAX(lambda, g%ub)
      end do
    end block
!
    if (PRESENT(edges)) then
      do concurrent(k=1:n_target - 1)
        call cantor_pair_inverse(core(k)%p, i, j)
        edges(1, k) = i
        edges(2, k) = j
      end do
    end if
    if (PRESENT(weights)) then
      do concurrent(k=1:n_target - 1)
        weights(k) = core(k)%lb
      end do
    end if
  end subroutine mobbrmsd_min_span_tree
!
! ---
!
! pure elemental function edge_to_ij(e) result(res)
!   type(edge), intent(in) :: e
!   type(ij)               :: res
!   call cantor_pair_inverse(e%p, res%i, res%j)
! end function edge_to_ij
!
  subroutine construct_chunk( &
 &             n_target, &
 &             n_core, &
 &             n_chunk, &
 &             header, &
 &             lambda, &
 &             par, &
 &             X, &
 &             n_pack, &
 &             e, &
 &             remove_com, &
 &             sort_by_g &
 &            )
    integer(IK), intent(in)        :: n_target
    integer(IK), intent(in)        :: n_core
    integer(IK), intent(in)        :: n_chunk
    integer(IK), intent(in)        :: par(n_target)
    type(mobbrmsd), intent(in)     :: header
    real(RK), intent(in)           :: lambda
    real(RK), intent(in)           :: X(*)
    integer(IK), intent(inout)     :: n_pack
    type(edge), intent(inout)      :: e(n_chunk)
    logical, intent(in), optional  :: remove_com, sort_by_g
    type(edge)                     :: t
    logical                        :: check1, check2
    integer(IK)                    :: i, j, k, l, m, c, n, ld, lc
!
    call init_edge(n_chunk, e)
!
    ld = n_target * (n_target - 1) / 2
    lc = MIN(n_chunk, ld - n_core)
!
    k = 0
    l = 0
    n_pack = 0
    do while (k < ld .and. l < lc)
      k = k + 1
      call cantor_pair_inverse(k, i, j)
      if (is_same(n_target, i, j, par)) cycle
      l = l + 1
      n_pack = n_pack + 1
    end do
!
    c = 0
    k = 0
    n = 0
!
    !$omp parallel shared(e,k,c), private(i,j,l,m,check1,check2)
    do
      !$omp critical
      l = k + 1
      m = c + 1
      if (k < ld .and. c < lc) then
        check1 = .false.
        call cantor_pair_inverse(l, i, j)
        k = l
        if (is_same(n_target, i, j, par)) then
          check2 = .true.
        else
          check2 = .false.
          c = m
          n = MAX(n, k)
        end if
      else
        check1 = .true.
      end if
      !$omp end critical
      if (check1) exit
      if (check2) cycle
      call update( &
   &             l, i, j, &
   &             header, &
   &             lambda, &
   &             X, &
   &             e(m), &
   &             remove_com, &
   &             sort_by_g &
   &            )
    end do
    !$omp end parallel
    do l = n_pack, 1, -1
      call down_heap(n_pack, l, e)
    end do
    k = n
    !$omp parallel shared(e,k), private(l,i,j,t,check1)
    do
      !$omp critical
      l = k + 1
      check1 = l > ld
      if (.not. check1) k = l
      !$omp end critical
      if (check1) exit
      call cantor_pair_inverse(l, i, j)
      if (is_same(n_target, i, j, par)) cycle
      call update( &
   &             l, i, j, &
   &             header, &
   &             lambda, &
   &             X, &
   &             t, &
   &             remove_com, &
   &             sort_by_g &
   &            )
      !$omp critical
      check1 = e(1)%lb < t%lb
      !$omp end critical
      if (check1) cycle
      !$omp critical
      e(1) = t
      call down_heap(n_chunk, 1, e)
      !$omp end critical
    end do
    !$omp end parallel
  end subroutine construct_chunk
!
  pure subroutine update_chunk( &
 &             n_target, &
 &             n_pack, &
 &             n_core, &
 &             n_chunk, &
 &             header, &
 &             lambda, &
 &             par, &
 &             X, &
 &             e, &
 &             remove_com, &
 &             sort_by_g &
 &            )
    integer(IK), intent(in)        :: n_target
    integer(IK), intent(in)        :: n_pack
    integer(IK), intent(in)        :: n_core
    integer(IK), intent(in)        :: n_chunk
    integer(IK), intent(in)        :: par(n_target)
    type(mobbrmsd), intent(in)     :: header
    real(RK), intent(in)           :: lambda
    real(RK), intent(in)           :: X(*)
    type(edge), intent(inout)      :: e(n_pack)
    logical, intent(in), optional  :: remove_com, sort_by_g
    real(RK)                       :: lambda_local
    type(edge)                     :: g
    integer(IK)                    :: i, j, k, p
    g = pick_heap(n_pack, n_target - n_core - 1, e)
    lambda_local = g%ub
    do while (lambda_local < lambda .and. .not. g%f)
      do k = 1, n_target - 1
        if (e(k)%f) cycle
        p = e(k)%p
        call cantor_pair_inverse(p, i, j)
        call update( &
     &             p, i, j, &
     &             header, &
     &             lambda_local, &
     &             X, &
     &             e(k), &
     &             remove_com, &
     &             sort_by_g &
     &            )
      end do
      do k = n_pack, 1, -1
        call down_heap(n_pack, k, e)
      end do
      g = pick_heap(n_pack, n_target - n_core - 1, e)
      lambda_local = MAX(lambda_local, g%ub)
    end do
  end subroutine update_chunk
!
  pure function pick_heap(n, p, h) result(res)
    integer(IK), intent(in) :: n, p
    type(edge), intent(in)  :: h(n)
    type(edge)              :: t(n), res
    integer(IK)             :: k
    t = h
    do k = n, n - p + 1, -1
      t(1) = h(k)
      call down_heap(k - 1, 1, t)
    end do
    res = t(1)
  end function pick_heap
!
  pure subroutine kruscal_core(n_target, n_pack, c, par, e, b)
    integer(IK), intent(in)    :: n_target, n_pack
    integer(IK), intent(inout) :: c, par(n_target)
    type(edge), intent(in)     :: e(n_pack)
    type(edge), intent(inout)  :: b(n_target - 1)
    type(edge)                 :: t(n_pack)
    integer(IK)                :: i, j, k
    logical                    :: is_same
    t = e
    do k = n_pack, 2, -1
      call swap(t(1), t(k))
      call down_heap(k - 1, 1, t)
    end do
    do k = 1, n_pack
      if (.not. t(k)%f) cycle
      call cantor_pair_inverse(t(k)%p, i, j)
      call unite(n_target, i, j, par, is_same)
      if (is_same) cycle
      c = c + 1
      b(c) = t(k)
      if (c >= n_target - 1) exit
    end do
  end subroutine kruscal_core
!
  pure elemental subroutine swap(a, b)
    type(edge), intent(inout) :: a, b
    type(edge)                :: s
    s = b
    b = a
    a = s
  end subroutine swap
!
  pure subroutine update( &
 &             k, i, j, &
 &             header, &
 &             lambda, &
 &             X, &
 &             e, &
 &             remove_com, &
 &             sort_by_g &
 &            )
    integer(IK), intent(in)       :: k, i, j
    type(mobbrmsd), intent(in)    :: header
    real(RK), intent(in)          :: lambda, X(*)
    type(edge), intent(inout)     :: e
    logical, intent(in), optional :: remove_com, sort_by_g
    type(mobbrmsd_state)          :: state
    real(RK)                      :: w(mobbrmsd_memsize(header))
    integer(IK)                   :: xpnt, ypnt, ldx
    ldx = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
    xpnt = (i - 1) * ldx + 1
    ypnt = (j - 1) * ldx + 1
    call mobbrmsd_run( &
   &       header, &
   &       state,  &
   &       X(xpnt), &
   &       X(ypnt), &
   &       W, &
   &       cutoff=lambda, &
   &       sort_by_g=sort_by_g &
   &      )
    e%w = 0
    e%f = mobbrmsd_state_is_finished(state)
    e%p = k
!   e%t = .false.
    e%lb = mobbrmsd_state_lowerbound_as_rmsd(state)
    e%ub = mobbrmsd_state_upperbound_as_rmsd(state)
  end subroutine update
!
  pure subroutine down_heap(n, p, e)
    integer(IK), intent(in)   :: n, p
    type(edge), intent(inout) :: e(*)
    integer(IK)               :: r, s, t
    r = p
    do while (r + r <= n)
      s = r + r
      t = s + 1
      if (t > n) then
        if (e(r)%lb > e(s)%lb) return
      else
        if (e(s)%lb > e(t)%lb) then
          if (e(r)%lb > e(s)%lb) return
        else
          if (e(r)%lb > e(t)%lb) return
          s = t
        end if
      end if
      call swap(e(r), e(s))
      r = s
    end do
  end subroutine down_heap
!
! pure recursive subroutine insert(p, q, a)
!   integer, intent(in)       :: p, q
!   type(edge), intent(inout) :: a(*)
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
!   type(edge), intent(inout)  :: e(*)
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
  pure function is_same(n, v1, v2, par)
    integer(IK), intent(in) :: n, v1, v2
    integer(IK), intent(in) :: par(n)
    logical                 :: is_same
    integer(IK)             :: r1, l1, r2, l2
    call findroot(n, v1, par, r1, l1)
    call findroot(n, v2, par, r2, l2)
    is_same = r1 == r2
  end function is_same
!
  pure subroutine unite(n, v1, v2, par, is_same)
    integer(IK), intent(in)    :: n, v1, v2
    integer(IK), intent(inout) :: par(n)
    logical, intent(inout)     :: is_same
    integer(IK)                :: r1, l1, r2, l2
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
  pure elemental subroutine cantor_pair_inverse(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    j = (INT(SQRT(real(8 * k - 7))) - 1) / 2 + 2 ! cantor_pair_inverse
    i = k - (j - 1) * (j - 2) / 2
  end subroutine cantor_pair_inverse
!
  pure subroutine init_edge(n, e)
    integer(IK), intent(in)   :: n
    type(edge), intent(inout) :: e(n)
    integer(IK)               :: i
    do concurrent(i=1:n)
      e(i)%p = 0
      e(i)%w = 0
      e(i)%f = .false.
!     e(i)%t = .false.
      e(i)%lb = HUGE(0.0_RK)
      e(i)%ub = HUGE(0.0_RK)
    end do
  end subroutine init_edge
!
end module mod_mobbrmsd_mst

