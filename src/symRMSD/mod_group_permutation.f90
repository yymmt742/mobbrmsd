!
!| Module for permutation.
module mod_group_permutation
  use mod_params, only: IK, RK
  implicit none
  private
  public :: group_permutation
  public :: group_permutation_tuple
  public :: group_permutation_swap
  public :: group_permutation_reverse
!
!| Module for permutation by decomposed cyclic groups. <br>
!  Swap indeces are stored in work array, in the following. <br>
!  - [[ p2, ..., pm ],  [ b1, b2, ..., bi, ..., bm ]] <br>
!  - pi (m-1)  :: pointer to bi. p1 is calculated by this%p + m - 1. <br>
!  - bi (m)    :: body. <br>
!  body is constructed as follows. <br>
!  - [n, s, s1, s2, ..., sj, ..., sn] <br>
!  - n  :: number of permutation. <br>
!  - s  :: order of cycle. <br>
!  - sj :: mapping sequence. <br>
!  - Total memsize is 2 + m*s
  type :: group_permutation
    private
    sequence
    !| pointers.
    integer(IK), public :: p = 1
    !| number of type cyclic orders.
    integer(IK)         :: m = 0
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
  pure function group_permutation_tuple_new(perm) result(res)
    integer(IK), intent(in)       :: perm(:)
    type(group_permutation_tuple) :: res
    integer(IK)                   :: n
!
    n = SIZE(perm)
!
!   early return
!
    if (n < 2) then
      res%t = group_permutation(1, 0)
      return
    elseif (n == 2) then
      if (ALL(perm(1:2) == [2, 1])) then
        res%t = group_permutation(1, 1)
        res%w = [1, 2, 2, 1]
      else
        res%t = group_permutation(1, 0)
      end if
      return
    elseif (n == 3) then
      if (ALL(perm(1:3) == [1, 3, 2])) then
        res%t = group_permutation(1, 1)
        res%w = [1, 2, 3, 2]
      elseif (ALL(perm(1:3) == [2, 1, 3])) then
        res%t = group_permutation(1, 1)
        res%w = [1, 2, 2, 1]
      elseif (ALL(perm(1:3) == [2, 3, 1])) then
        res%t = group_permutation(1, 1)
        res%w = [1, 3, 2, 3, 1]
      elseif (ALL(perm(1:3) == [3, 1, 2])) then
        res%t = group_permutation(1, 1)
        res%w = [1, 3, 3, 1, 2]
      elseif (ALL(perm(1:3) == [3, 2, 1])) then
        res%t = group_permutation(1, 1)
        res%w = [1, 2, 3, 1]
      else
        res%t = group_permutation(1, 0)
      end if
      return
    elseif (is_not_permutation(n, perm)) then
      res%t = group_permutation(1, 0)
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
      res%t = group_permutation(1, k)
      if (res%t%m < 1) return
      j = res%t%m - 1
      do i = 1, k
        j = j + w(n + l + i) * w(n + l + k + i) + 2
      end do
      ALLOCATE(res%w(j))
!
      res%w(k) = w(n + l + k + 1)
      res%w(k + 1) = w(n + l + 1)
      call proc(n, l, w(n + l + 1), w, w(n + 1), res%w(k+2))
      j = k
      do i = 2, k
        j = j + 2 + w(n + l + i - 1) * w(n + l + k + i - 1)
        res%w(i - 1) = j
        res%w(j) = w(n + l + k + i)
        res%w(j + 1) = w(n + l + i)
        call proc(n, l, w(n + l + i), w, w(n + 1), res%w(j+2))
      end do
!
    end block
!
  end function group_permutation_tuple_new
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
! w is now stored [b1, b2, ..., bl, s1, s2, ..., sl, t1, ..., tk, n1, ..., nk]
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
  pure subroutine group_permutation_swap(this, w, d, X)
    type(group_permutation), intent(in) :: this
    !! group_permutation
    integer(IK), intent(in)             :: w(*)
    !! work array.
    integer(IK), intent(in)             :: d
    !! leading dimension of x
    real(RK), intent(inout)             :: X(*)
    !! data array.
    integer(IK)                         :: i
!
    do concurrent(i=1:this%m)
      block
        integer(IK) :: j, p, n, s, ns
        if (i == 1) then
          p = this%p + this%m
        else
          p = this%p + w(this%p + i - 2)
        end if
        n = w(p - 1)
        s = w(p)
        ns = p + n * s
        do concurrent(j=p + 1:ns:s)
          call cyclic_swap_real(d, s, w(j), X)
        end do
      end block
    end do
!
  end subroutine group_permutation_swap
!
  pure subroutine group_permutation_reverse(this, w, d, X)
    type(group_permutation), intent(in) :: this
    !! group_permutation
    integer(IK), intent(in)             :: w(*)
    !! work array.
    integer(IK), intent(in)             :: d
    !! leading dimension of x
    real(RK), intent(inout)             :: X(*)
    !! data array.
    integer(IK)                         :: i
!
    do concurrent(i=1:this%m)
      block
        integer(IK) :: j, p, n, s, ns
        if (i == 1) then
          p = this%p + this%m
        else
          p = this%p + w(this%p + i - 2)
        end if
        n = w(p - 1)
        s = w(p)
        ns = p + n * s
        do concurrent(j=p + 1:ns:s)
          call cyclic_reverse_real(d, s, w(j), X)
        end do
      end block
    end do
!
  end subroutine group_permutation_reverse
!
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
  pure subroutine cyclic_reverse_real(d, s, q, X)
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
  end subroutine cyclic_reverse_real
!
  pure elemental subroutine group_permutation_tuple_destroy(this)
    type(group_permutation_tuple), intent(inout) :: this
    if(ALLOCATED(this%w)) deallocate(this%w)
  end subroutine group_permutation_tuple_destroy
!
end module mod_group_permutation

