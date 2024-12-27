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
    !  integer(IK) :: l, r, p, s, q, w
    logical     :: f, t
    real(RK)    :: lb
  end type
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
    integer(IK)                         :: f(n_target), pi, pj
    type(edge)                          :: e(n_target * (n_target - 1) / 2)
    type(edge)                          :: b(n_target - 1)
    real(RK)                            :: ub, rmsd, cutoff_global
    !integer(IK)                         :: root, leaf
    integer(IK)                         :: n_chunk
    integer(IK)                         :: i, j, k, c, ilcl, jlcl, nlim, spnt, xpnt, ypnt, wpnt, l0, ldx, ldw
!
    !n_chunk = n_target * (n_target - 1) / 2
    n_chunk = MIN(n_target * 3, n_target * (n_target - 1) / 2)
    cutoff_global = 10.0_RK
    do
      print *, cutoff_global
      call construct_chunk( &
         &  n_target, &
         &  n_chunk, &
         &  header, &
         &  cutoff_global, &
         &  X, &
         &  e, &
         &  remove_com, &
         &  sort_by_g &
         & )
      print'(10f9.3)', e(:n_chunk)%lb
!     call kruscal( &
!        &  n_target, &
!        &  header, &
!        &  cutoff_global, &
!        &  X, &
!        &  e, &
!        &  root, &
!        &  leaf, &
!        &  remove_com, &
!        &  sort_by_g &
!        & )
!     print *, root, leaf
!     k = root
!     do
!       call cantor_pair_inverse(k, i, j)
!       print'(3I8,L4,2f6.3)', i, j, k, e(k)%f, e(k)%lb, e(k)%ub
!       k = e(k)%q
!       if (k < 1) exit
!     end do
!     print *
!     if (e(leaf)%f) exit
!     cutoff_global = e(leaf)%ub * 1.0001
      call kruscal(n_target, n_chunk, c, e, b)
      exit
    end do
!
    if (PRESENT(edges)) then
      do concurrent(k=1:n_target - 1)
        call cantor_pair_inverse(b(k)%p, i, j)
        edges(1, k) = i
        edges(2, k) = j
      end do
    end if
    if (PRESENT(weights)) then
      do concurrent(k=1:n_target - 1)
        weights(k) = b(k)%lb
      end do
    end if
    return
!
    l0 = n_target * (n_target - 1) / 2
    ldx = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
    ldw = mobbrmsd_memsize(header)
!
    !$omp parallel do private(spnt, xpnt, ypnt, wpnt)
    do j = 1, n_target
      do i = 1, j - 1
        call cantor_pair(i, j, spnt)
        xpnt = (i - 1) * ldx + 1
        ypnt = (j - 1) * ldx + 1
        wpnt = (spnt - 1) * ldw + 1
        call mobbrmsd_run( &
       &       header, &
       &       state(spnt), &
       &       X(xpnt), &
       &       X(ypnt), &
       &       W(wpnt), &
       &       cutoff=cutoff, &
       &       ub_cutoff=ub_cutoff, &
       &       difflim=difflim, &
       &       maxeval=1, &
       &       remove_com=remove_com, &
       &       sort_by_g=sort_by_g, &
       &       difflim_absolute=difflim_absolute  &
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
    do concurrent(i=1:n_target)
      f(i) = i
    end do
    k = MINLOC(mobbrmsd_state_lowerbound_as_rmsd(state(:l0)), 1)
    call cantor_pair_inverse(k, i, j)
    pi = 0
    pj = 0
    f(1) = f(i)
    f(i) = 1
!
    do k = 1, n_target - 1
      ub = cutoff_global
      do i = 1, k
        do j = k + 1, n_target
          call cantor_pair(f(i), f(j), spnt)
          ub = MIN(mobbrmsd_state_rmsd(state(spnt)), ub)
        end do
      end do
      nlim = k * (n_target - k)
      i = 0
      !$omp parallel private(ilcl, jlcl, spnt, xpnt, ypnt, wpnt)
      do
        !$omp critical
        ilcl = i
        i = i + 1
        !$omp end critical
        if (ilcl >= nlim) exit
        jlcl = 1 + MODULO(ilcl, k)
        ilcl = k + 1 + ilcl / k
        call cantor_pair(f(ilcl), f(jlcl), spnt)
        if (.not. mobbrmsd_state_is_finished(state(spnt))) then
          wpnt = (spnt - 1) * ldw + 1
          call mobbrmsd_restart( &
         &       header, &
         &       state(spnt), &
         &       W(wpnt), &
         &       cutoff=ub, &
         &       ub_cutoff=ub_cutoff, &
         &       difflim=difflim, &
         &       maxeval=maxeval &
         &      )
        end if
        !$omp critical
        ub = MIN(ub, mobbrmsd_state_rmsd(state(spnt)))
        !$omp end critical
      end do
      !$omp end parallel
!
      ub = RHUGE
      do i = 1, k
        do j = k + 1, n_target
          call cantor_pair(f(i), f(j), spnt)
          rmsd = mobbrmsd_state_rmsd(state(spnt))
          if (rmsd < ub) then
            ub = rmsd
            pi = i
            pj = j
          end if
        end do
      end do
!
      if (PRESENT(edges)) then
        if (f(pi) < f(pj)) then
          edges(1, k) = f(pi)
          edges(2, k) = f(pj)
        else
          edges(1, k) = f(pj)
          edges(2, k) = f(pi)
        end if
      end if
      if (PRESENT(weights)) then
        weights(k) = ub
      end if
      j = f(k + 1)
      f(k + 1) = f(pj)
      f(pj) = j
    end do
!
  end subroutine mobbrmsd_min_span_tree
!
! ---
!
  subroutine construct_chunk( &
 &             n_target, &
 &             n_chunk, &
 &             header, &
 &             cutoff, &
 &             X, &
 &             e, &
 &             remove_com, &
 &             sort_by_g &
 &            )
    integer(IK), intent(in)        :: n_target
    integer(IK), intent(in)        :: n_chunk
    type(mobbrmsd), intent(in)     :: header
    real(RK), intent(in)           :: cutoff
    real(RK), intent(in)           :: X(*)
    type(edge), intent(inout)      :: e(n_chunk)
    logical, intent(in), optional  :: remove_com
    logical, intent(in), optional  :: sort_by_g
    type(edge)                     :: t
    integer(IK)                    :: k
    !$omp parallel do
    do k = 1, n_chunk
      call update( &
   &             k, &
   &             header, &
   &             cutoff, &
   &             X, &
   &             e(k), &
   &             remove_com, &
   &             sort_by_g &
   &            )
    end do
    !$omp end parallel do
    do k = n_chunk, 1, -1
      call down_heap(n_chunk, k, e)
    end do
    !$omp parallel do shared(e), private(k,t)
    do k = n_chunk + 1, n_target * (n_target - 1) / 2
      call update( &
   &             k, &
   &             header, &
   &             cutoff, &
   &             X, &
   &             t, &
   &             remove_com, &
   &             sort_by_g &
   &            )
      if (e(1)%lb < t%lb) cycle
      !$omp critical
      e(1) = t
      call down_heap(n_chunk, 1, e)
      !$omp end critical
    end do
    !$omp end parallel do
    do k = n_chunk, 2, -1
      call swap(e(1), e(k))
      call down_heap(k - 1, 1, e)
    end do
  end subroutine construct_chunk
!
  pure subroutine kruscal(n_target, n_chunk, c, e, b)
    integer(IK), intent(in)    :: n_target, n_chunk
    integer(IK), intent(inout) :: c
    type(edge), intent(inout)  :: e(n_chunk), b(n_target - 1)
    integer(IK)                :: par(n_target)
    integer(IK)                :: i, j, k
    call init_par(n_target, par)
    c = 0
    do k = 1, n_chunk
      call cantor_pair_inverse(e(k)%p, i, j)
      call unite(n_target, i, j, par, e(k)%t)
      if (e(k)%t) cycle
      c = c + 1
      b(c) = e(k)
      if (c >= n_target - 1) exit
    end do
  end subroutine kruscal
!
! subroutine kruscal( &
!&             n_target, &
!&             header, &
!&             cutoff, &
!&             X, &
!&             e, &
!&             root, &
!&             leaf, &
!&             remove_com, &
!&             sort_by_g &
!&            )
!   integer(IK), intent(in)        :: n_target
!   type(mobbrmsd), intent(in)     :: header
!   real(RK), intent(in)           :: cutoff
!   real(RK), intent(in)           :: X(*)
!   type(edge), intent(inout)      :: e(n_target * (n_target - 1) / 2)
!   integer(IK), intent(inout)     :: root
!   integer(IK), intent(inout)     :: leaf
!   logical, intent(in), optional  :: remove_com
!   logical, intent(in), optional  :: sort_by_g
!   integer(IK)                    :: heap(n_target * (n_target - 1) / 2)
!   integer(IK)                    :: par(n_target), i, j, k, l, c
!   call init_par(SIZE(heap), heap)
!   do k = SIZE(e), 1, -1
!     call update( &
!  &             k, &
!  &             header, &
!  &             cutoff, &
!  &             X, &
!  &             e(k), &
!  &             remove_com, &
!  &             sort_by_g &
!  &            )
!     call down_heap(SIZE(e), k, e, heap)
!   end do
!   do k = SIZE(heap), 2, -1
!     call swap(heap(1), heap(k))
!     call down_heap(k - 1, 1, e, heap)
!   end do
!   root = heap(1)
!   do concurrent(k=2:SIZE(e))
!     e(heap(k - 1))%s = heap(k)
!   end do
!   call init_par(n_target, par)
!   k = root
!   l = root
!   c = 0
!   do while (c < n_target - 1)
!     call cantor_pair_inverse(k, i, j)
!     call unite(n_target, i, j, par, e(k)%t)
!     if (e(k)%t) then
!       k = e(k)%s
!       cycle
!     end if
!     e(l)%q = k
!     l = k
!     k = e(k)%s
!     c = c + 1
!   end do
!   leaf = l
!
! end subroutine kruscal
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
 &             k, &
 &             header, &
 &             cutoff, &
 &             X, &
 &             e, &
 &             remove_com, &
 &             sort_by_g &
 &            )
    integer(IK), intent(in)             :: k
    type(mobbrmsd), intent(in)          :: header
    real(RK), intent(in)                :: cutoff
    real(RK), intent(in)                :: X(*)
    type(edge), intent(inout)           :: e
    logical, intent(in), optional       :: remove_com
    logical, intent(in), optional       :: sort_by_g
    type(mobbrmsd_state)                :: state
    real(RK)                            :: w(mobbrmsd_memsize(header))
    integer                             :: xpnt, ypnt
    integer                             :: i, j, ldx
    ldx = mobbrmsd_n_dims(header) * mobbrmsd_n_atoms(header)
    call cantor_pair_inverse(k, i, j)
    xpnt = (i - 1) * ldx + 1
    ypnt = (j - 1) * ldx + 1
    call mobbrmsd_run( &
   &       header, &
   &       state,  &
   &       X(xpnt), &
   &       X(ypnt), &
   &       W, &
   &       cutoff=cutoff, &
   &       sort_by_g=sort_by_g &
   &      )
    !e%l = 0
    !e%r = 0
    !e%s = 0
    !e%q = 0
    e%w = 0
    e%p = k
    e%f = mobbrmsd_state_is_finished(state)
    e%t = .false.
    e%lb = mobbrmsd_state_lowerbound_as_rmsd(state)
    !e%ub = mobbrmsd_state_upperbound_as_rmsd(state)
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
  pure elemental subroutine cantor_pair_inverse(k, i, j)
    integer(IK), intent(in)    :: k
    integer(IK), intent(inout) :: i, j
    j = (INT(SQRT(real(8 * k - 7))) - 1) / 2 + 2 ! cantor_pair_inverse
    i = k - (j - 1) * (j - 2) / 2
  end subroutine cantor_pair_inverse
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
  pure elemental subroutine cantor_pair(i, j, k)
    integer(IK), intent(in)    :: i, j
    integer(IK), intent(inout) :: k
    if (i > j) then
      k = (i - 2) * (i - 1) / 2 + j
    else
      k = (j - 2) * (j - 1) / 2 + i
    end if
  end subroutine cantor_pair
end module mod_mobbrmsd_mst

