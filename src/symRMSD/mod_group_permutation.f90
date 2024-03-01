!
!| Module for permutation.
module mod_group_permutation
  use mod_params, only: IK, RK
  implicit none
  private
  public :: group_permutation
  public :: group_permutation_tuple
  public :: group_permutation_swap
  public :: group_permutation_inverse
!
!| Module for permutation by decomposed cyclic groups. <br>
!  gp are stored in work array, in the following. <br>
!  - [p2, ..., ps, swp_1, swp_2, ..., swp_s]
!  - pi (s-1)  :: pointer to swp_2 to swp_s. <br>
!  - p1 is calculated in p + s - 1 <br>
!  - pi is calculated in p + w(i-2) <br>
!  Swap indeces are stored in work array, in the following. <br>
!  - [[p1, p2, ..., pm],  [ b1, b2, ..., bi, ..., bm ]] <br>
!  - m         :: number of order - 1. <br>
!  - pi (m)    :: pointer to bi+1. Note that m = p1 - 1. <br>
!  - bi (m+1)  :: body. <br>
!  body is constructed as follows. <br>
!  - [n, s, s1, s2, ..., sj, ..., sn] <br>
!  - n  :: number of permutation. <br>
!  - s  :: order of cycle. <br>
!  - sj :: mapping sequence. <br>
!  - Total memsize is 2 + m * s
  type :: group_permutation
    private
    sequence
    !| pointers.
    integer(IK), public :: p = 1
    !| number of permutation. if s==0, gp only has identity mapping.
    integer(IK)         :: s = 0
  end type group_permutation
!
!| A set of t and w arrays. <br>
!  This is mainly used for passing during initialization.
  type :: group_permutation_tuple
    !| group_permutation.
    type(group_permutation)  :: t
    !| w array.
    integer(IK), allocatable :: w(:)
  contains
    final :: group_permutation_tuple_destroy
  end type group_permutation_tuple
!
  interface group_permutation_tuple
    module procedure group_permutation_tuple_new
  end interface
!
contains
!
!| Constructor. <br>
  function group_permutation_tuple_new(perm) result(res)
    integer(IK), intent(in)       :: perm(:, :)
    !! codomains, [[a1,a2,...,am],[b1,b2,...,bm],...].
    type(group_permutation_tuple) :: res
    !! return value.
    integer(IK), allocatable      :: t(:)
    integer(IK)                   :: i, n, s, p
!
    n = SIZE(perm, 1)
    s = SIZE(perm, 2)
!
    res%t = group_permutation(1, s)
    allocate (res%w(MAX(s - 1, 0)))
!
    if (s < 1) return
!
    p = 0
    t = decompose_to_cyclic(perm(:, 1))
    res%w = [res%w, t]
!
    do i = 2, s
      p = p + SIZE(t)
      res%w(i - s) = p
      t = decompose_to_cyclic(perm(:, 1))
      res%w = [res%w, t]
    end do
!
  end function group_permutation_tuple_new
!
  function decompose_to_cyclic(perm) result(res)
    integer(IK), intent(in)  :: perm(:)
    integer(IK), allocatable :: res(:)
    integer(IK)              :: n
!
!   early return
!
    n = SIZE(perm)
    ALLOCATE(res, source=[0])
!
    if (n < 2) then
      return
    elseif (n == 2) then
      if (ALL(perm(1:2) == [2, 1])) res = [1, 1, 2, 2, 1]
      return
    elseif (n == 3) then
      if (ALL(perm(1:3) == [1, 3, 2])) then
        res = [1, 1, 2, 3, 2]
      elseif (ALL(perm(1:3) == [2, 1, 3])) then
        res = [1, 1, 2, 2, 1]
      elseif (ALL(perm(1:3) == [2, 3, 1])) then
        res = [1, 1, 3, 2, 3, 1]
      elseif (ALL(perm(1:3) == [3, 1, 2])) then
        res = [1, 1, 3, 3, 1, 2]
      elseif (ALL(perm(1:3) == [3, 2, 1])) then
        res = [1, 1, 2, 3, 1]
      end if
      return
    elseif (is_not_permutation(n, perm)) then
      return
    end if
!
!   decompose to cyclic permutations.
!
    block
      integer(IK)                   :: i, j, k, l, w(4 * n)
!
      do concurrent(i=1:n + n)
        w(i) = -1
      end do
      do concurrent(i=1:n)
        w(n + n + i) = perm(i)
      end do
!
      i = 1
      j = n + 1
      k = n + n
!
      do l = 1, n
        k = k + 1
        if (w(k) < 0) cycle
        call count_permutation(n, w(n + n + 1), l, w(i), w(j))
        i = i + w(j)
        j = j + 1
      end do
!
      l = j - n - 1
      k = n + l
!
!     w(1:n) is now stored [b1, b2, ..., bl]
!     w(n+1:n+l) is now stored [s1, s2, ..., sl]
!
      do concurrent(i=1:l)
        w(k + i) = w(n + i)
      end do
      call uniq(l, k, w(n + l + 1))
!
!     k                  :: n, number of uniq cyclic order (without 1).
!     w(n+l+1:n+l+k)     :: [t1, ..., tk], uniq cyclic order list.
!     w(n+l+k+1:n+l+k+k) :: [n1, ..., nk], count of ti.
!
!     w is now stored [b1, b2, ..., bl, s1, s2, ..., sl, t1, ..., tk, n1, ..., nk]
!
      if (k < 1) return
!
!     count array size
      j = 3 * k + SUM([(w(n + l + i) * w(n + l + k + i), i=1,k)])
!
!     pack w to array.
      res = proc_w(n, j, k, l, w(1), w(n + 1), w(n + l + k + 1), w(n + l + 1))
!
    end block
!
  end function decompose_to_cyclic
!
  pure function proc_w(nd, nj, nk, nl, b, s, n, t) result(res)
    integer(IK), intent(in) :: nd, nj, nk, nl
    integer(IK), intent(in) :: n(nk), t(nk), b(nd), s(nl)
    integer(IK)             :: res(nj)
    integer(IK)             :: i, j
!
      res(1) = nk + 1
      res(nk + 1) = n(1)
      res(nk + 2) = t(1)
      call proc(nd, nl, t(1), b, s(1), res(nk + 3))
      j = nk + 1
      do i = 2, nk
        j = j + 2 + t(i - 1) * n(i - 1)
        res(i) = j
        res(j) = n(i)
        res(j + 1) = t(i)
        call proc(nd, nl, t(i), b, s, res(j + 2))
      end do
!
  end function proc_w
!
  pure subroutine proc(nd, l, t, b, s, res)
    integer(IK), intent(in)    :: nd, l, t, b(nd), s(l)
    integer(IK), intent(inout) :: res(*)
    integer(IK)                :: i, j, k
!
    j = 1
    k = 1
!
    do i = 1, l
      if (t==s(i))then
        res(j:j + t - 1) = b(k:k + t - 1)
        j = j + s(i)
      endif
      k = k + s(i)
    end do
!
  end subroutine proc
!
  pure function is_not_permutation(n, perm) result(res)
    integer(IK), intent(in) :: n, perm(*)
    logical                 :: res
    integer(IK)             :: i
    do i = 1, n - 1
      res = perm(i) < 1 .or. n < perm(i) .or. ANY(perm(i) == perm(i + 1:n))
      if(res) return
    end do
    res = perm(n) < 1 .or. n < perm(n)
  end function is_not_permutation
!
!| returns cyclic indices and its order.
  pure subroutine count_permutation(n, perm, s, l, c)
    integer(IK), intent(in)    :: n, s
    integer(IK), intent(inout) :: perm(*), l(*), c
    integer(IK)                :: e0, e
!
    c = 1
    l(1) = s
    e0 = s
    e = perm(s)
    perm(s) = -perm(s)
!
    do
      if (e < 1 .or. n < e) then
        c = -1; return
      end if
      if (e0 == e) return
      c = c + 1
      l(c) = e
      e = perm(e)
      perm(l(c)) = - perm(l(c))
    end do
!
  end subroutine count_permutation
!
!| returns uniq integer list without 1.
  pure subroutine uniq(n, u, res)
    integer(IK), intent(in)    :: n
    integer(IK), intent(inout) :: u, res(*)
    integer(IK)                :: i, j
!
    u = 0
    if (n > 1) call qsi(n, res)
    if(res(1) > 1) u = 1
!
    do concurrent(i=n + 1 + u:n + n + 1)
      res(i) = 0
    end do
!
    j = n + 1 + u
!
    do i = 1, n - 1
      res(j) = res(j) + 1
      if (res(i) == res(i + 1)) cycle
      u = u + 1
      j = j + 1
      res(u) = res(i + 1)
    end do
!
    res(j) = res(j) + 1
!
    do concurrent(i=1:u)
      res(u + i) = res(n + i + 1)
    end do
!
  end subroutine uniq
!
  pure recursive subroutine qsi(n, s)
    integer(IK), intent(in)    :: n
    integer(IK), intent(inout) :: s(n)
    integer(IK)                :: i, j, p, t
!
    j = n; i = 1; p = j / 2
!
    do
!
      do while (s(i) < s(p)); i = i + 1; end do
      do while (s(j) > s(p)); j = j - 1; end do
!
      if (i >= j) exit
!
      t = s(i); s(i) = s(j); s(j) = t
      if (i == p) then; p = j
      elseif (j == p) then; p = i
      end if
      i = i + 1; j = j - 1
!
    end do
!
    if (2 < i) call qsi(i - 1, s(:i - 1))
    if (j + 1 < n) call qsi(n - j, s(j + 1:))
!
  end subroutine qsi
!
! pure subroutine group_permutation_swap_int_l1(this, X)
!   class(group_permutation), intent(inout) :: this
!   integer(IK), intent(inout)              :: X(*)
!   integer(IK)                             :: i, j
!
!   if (.not. ALLOCATED(this%p)) return
!
!   do concurrent(i=1:SIZE(this%p, 2))
!     do concurrent(j=1:this%p(3, i))
!       block
!         integer(IK) :: b
!         b = this%p(1, i) + this%p(2, i) * (j - 1)
!         call cyclic_swap_int(1, this%p(2, i), this%q(b), X)
!       end block
!     end do
!   end do
!
! end subroutine group_permutation_swap_int_l1
!
! pure subroutine group_permutation_swap_int_l2(this, d, X)
!   class(group_permutation), intent(inout) :: this
!   integer(IK), intent(in)                 :: d
!   integer(IK), intent(inout)              :: X(d, *)
!   integer(IK)                             :: i, j
!
!   if (.not. ALLOCATED(this%p)) return
!
!   do concurrent(i=1:SIZE(this%p, 2))
!     do concurrent(j=1:this%p(3, i))
!       block
!         integer(IK) :: b
!         b = this%p(1, i) + this%p(2, i) * (j - 1)
!         call cyclic_swap_int(d, this%p(2, i), this%q(b), X)
!       end block
!     end do
!   end do
!
! end subroutine group_permutation_swap_int_l2
!
   subroutine group_permutation_swap(this, w, s, d, X)
    type(group_permutation), intent(in) :: this
    !! group_permutation
    integer(IK), intent(in) :: w(*)
    !! work array.
    integer(IK), intent(in) :: s
    !! permutation index.
    integer(IK), intent(in) :: d
    !! leading dimension of x
    real(RK), intent(inout) :: X(*)
    !! data array.
!
    if (s < 1 .or. this%s < s) return ! identity map
    call swap_real(w(this%p + s - 1), d, X)
!
  end subroutine group_permutation_swap
!
  pure subroutine swap_real(w, d, X)
    integer(IK), intent(in) :: w(*)
    !! swap indices.
    integer(IK), intent(in) :: d
    !! leading dimension of x
    real(RK), intent(inout) :: X(*)
    !! data array.
    integer(IK)             :: i, m
!
    m = w(1) - 1
!
    do concurrent(i=1:m)
      block
        integer(IK) :: j, p, n, s, ns
        p = w(i)
        n = w(p)
        p = p + 1
        s = w(p)
        p = p + 1
        ns = p + (n - 1) * s
        do concurrent(j=p:ns:s)
          call cyclic_swap_real(d, s, w(j), X)
        end do
      end block
    end do
!
  end subroutine swap_real
!
  pure subroutine group_permutation_inverse(this, w, s, d, X)
    type(group_permutation), intent(in) :: this
    !! group_permutation
    integer(IK), intent(in) :: w(*)
    !! work array.
    integer(IK), intent(in) :: s
    !! permutation index.
    integer(IK), intent(in) :: d
    !! leading dimension of x
    real(RK), intent(inout) :: X(*)
    !! data array.
!
    if (s < 1 .or. this%s < s) return ! identity map
    call inverse_real(w(this%p + s - 1), d, X)
!
  end subroutine group_permutation_inverse
!
  pure subroutine inverse_real(w, d, X)
    integer(IK), intent(in) :: w(*)
    !! swap indices.
    integer(IK), intent(in) :: d
    !! leading dimension of x
    real(RK), intent(inout) :: X(*)
    !! data array.
    integer(IK)             :: i, m
!
    m = w(1) - 1
!
    do concurrent(i=1:m)
      block
        integer(IK) :: j, p, n, s, ns
        p = w(i)
        n = w(p)
        p = p + 1
        s = w(p)
        p = p + 1
        ns = p + (n - 1) * s
        do concurrent(j=p:ns:s)
          call cyclic_inverse_real(d, s, w(j), X)
        end do
      end block
    end do
!
  end subroutine inverse_real
!
! pure subroutine cyclic_swap_int(d, s, q, X)
!   integer(IK), intent(in)    :: d, s, q(*)
!   integer(IK), intent(inout) :: X(d, *)
!   integer(IK)                :: T(d)
!   integer(IK)                :: i, j
!   do concurrent(i=1:d)
!     T(i) = X(i, q(1))
!   end do
!   do j = 1, s - 1
!     do concurrent(i=1:d)
!       X(i, q(j)) = X(i, q(j + 1))
!     end do
!   end do
!   do concurrent(i=1:d)
!     X(i, q(s)) = T(i)
!   end do
! end subroutine cyclic_swap_int
!
  pure subroutine cyclic_swap_real(d, s, q, X)
    integer(IK), intent(in) :: d, s, q(*)
    real(RK), intent(inout) :: X(d, *)
    real(RK)                :: T(d)
    integer(IK)             :: i, j
    do concurrent(i=1:d)
      T(i) = X(i, q(1))
    end do
    do j = 1, s - 1
      do concurrent(i=1:d)
        X(i, q(j)) = X(i, q(j + 1))
      end do
    end do
    do concurrent(i=1:d)
      X(i, q(s)) = T(i)
    end do
  end subroutine cyclic_swap_real
!
  pure subroutine cyclic_inverse_real(d, s, q, X)
    integer(IK), intent(in) :: d, s, q(*)
    real(RK), intent(inout) :: X(d, *)
    real(RK)                :: T(d)
    integer(IK)             :: i, j
    do concurrent(i=1:d)
      T(i) = X(i, q(s))
    end do
    do j = s - 1, 1, -1
      do concurrent(i=1:d)
        X(i, q(j + 1)) = X(i, q(j))
      end do
    end do
    do concurrent(i=1:d)
      X(i, q(1)) = T(i)
    end do
  end subroutine cyclic_inverse_real
!
  pure elemental subroutine group_permutation_tuple_destroy(this)
    type(group_permutation_tuple), intent(inout) :: this
    if(ALLOCATED(this%w)) deallocate(this%w)
  end subroutine group_permutation_tuple_destroy
!
end module mod_group_permutation

