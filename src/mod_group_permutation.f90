module mod_group_permutation
  use mod_params, only: IK, RK
  implicit none
  private
  public :: group_permutation
!
  integer(IK), parameter :: DEF_n = 1
  !! default number of components.
!
  type :: group_permutation
    private
    integer(IK), allocatable :: p(:, :)
    integer(IK), allocatable :: q(:)
  contains
    procedure          :: nfree        => group_permutation_nfree
    procedure          :: free_indices => group_permutation_free_indices
    procedure, private :: group_permutation_swap_int_l1
    procedure, private :: group_permutation_swap_int_l2
    procedure, private :: group_permutation_swap_real_l1
    procedure, private :: group_permutation_swap_real_l2
    generic            :: swap         => group_permutation_swap_int_l1, &
                      &                   group_permutation_swap_int_l2, &
                      &                   group_permutation_swap_real_l1, &
                      &                   group_permutation_swap_real_l2
    procedure          :: clear        => group_permutation_clear
    final              :: group_permutation_destroy
  end type group_permutation
!
  interface group_permutation
    module procedure group_permutation_new
  end interface group_permutation
!
contains
!
!| Constructer
  pure function group_permutation_new(prm) result(res)
    integer(IK), intent(in) :: prm(:)
    !! permutation indices
    type(group_permutation) :: res
!
    allocate (res%p(0, 0))
    allocate (res%q(0))
    call decompose_to_cyclic(SIZE(prm), prm, res%p, res%q)
!
  end function group_permutation_new
!
  pure subroutine decompose_to_cyclic(n, prm, p, q)
    integer(IK), intent(in)                 :: n, prm(n)
    integer(IK), intent(inout), allocatable :: p(:, :), q(:)
    integer(IK)                             :: w(5 * n)
    integer(IK)                             :: i, j, k, l
!
    if (n < 1) return
!
    w(1) = -1
    w(:3*n) = -1
    do k = 1, n
      if (prm(k) < 1 .or. n < prm(k)) return
      w(k) = prm(k)
    end do
!
    i = n + n + 1
    j = n + 1
    k = 0
!
    do l = 1, n
      if (w(l) < 0) cycle
      call count_permutation(n, w, l, w(i), w(j))
      if (w(j) < 0) then
        w(1) = -1; return
      end if
      i = i + w(j)
      j = j + 1
    enddo
!
    w(1) = j - n - 1
    do concurrent(i=1:w(1))
      w(i + 1) = w(i + n)
    end do
    do concurrent(i=1:n)
      w(i + w(1) + 1) = w(n + n + i)
    end do
!
    i = w(1) + n + 2
    j = w(1) + 2
    call uniq(w(1), w(2), w(i), w(i + 1))
    k = i + w(i) + w(i) + 1
    call proc(n, w(1), w(i), w(2), w(j), w(i + 1), w(i + w(i) + 1), w(k), w(1))
!
    p = header(w(1), w(2))
    i = 2 * w(1) + 2
    j = i + p(1, w(1)) + p(2, w(1)) * p(3, w(1)) - 2
    q = w(i:j)
!
  end subroutine decompose_to_cyclic
!
  pure function header(s, w) result(res)
    integer(IK), intent(in)    :: s, w(2, *)
    integer(IK)                :: res(3, s)
    integer(IK)                :: i
    do concurrent(i=1:s)
      res(2, i) = w(1, i)
      res(3, i) = w(2, i)
    end do
    res(1, 1) = 1
    do i = 1, s - 1
      res(1, i + 1) = res(1, i) + res(2, i) * res(3, i)
    end do
  end function header
!
  pure subroutine proc(n, m, s, p, q, u, t, w, res)
    integer(IK), intent(in)    :: n, m, s, p(*), q(*), u(*), t(*)
    integer(IK), intent(inout) :: w(*), res(*)
    integer(IK)                :: i, j, k, l
!
    do concurrent(i=1:s)
      w(2 * i - 1) = u(i)
      w(2 * i)     = t(i)
    end do
!
    l = s + s + 1
!
    do i = 1, s
      k = 0
      do j = 1, m
        k = k + p(j)
        if (p(j) /= u(i)) cycle
        w(l:l + u(i) - 1) = q(k - u(i) + 1:k)
        l = l + u(i)
      end do
    end do
!
    res(1) = s
    do concurrent(i=1:s + s + n)
      res(i + 1) = w(i)
    end do
!
  end subroutine proc
!
  pure subroutine count_permutation(n, prm, s, l, c)
    integer(IK), intent(in)    :: n, s
    integer(IK), intent(inout) :: prm(*), l(*), c
    integer(IK)                :: e0, e
!
    c = 1
    l(1) = s
    e0 = s
    e = prm(s)
    prm(s) = -prm(s)
!
    do
      if (e < 1 .or. n < e) then
        c = -1; return
      end if
      if (e0 == e) return
      c = c + 1
      l(c) = e
      e = prm(e)
      prm(l(c)) = - prm(l(c))
    end do
!
  end subroutine count_permutation
!
  pure subroutine group_permutation_swap_int_l1(this, X)
    class(group_permutation), intent(inout) :: this
    integer(IK), intent(inout)              :: X(*)
    integer(IK)                             :: i, j
!
    if (.not. ALLOCATED(this%p)) return
!
    do concurrent(i=1:SIZE(this%p, 2))
      do concurrent(j=1:this%p(3, i))
        block
          integer(IK) :: b
          b = this%p(1, i) + this%p(2, i) * (j - 1)
          call cyclic_swap_int(1, this%p(2, i), this%q(b), X)
        end block
      end do
    end do
!
  end subroutine group_permutation_swap_int_l1
!
  pure subroutine group_permutation_swap_int_l2(this, d, X)
    class(group_permutation), intent(inout) :: this
    integer(IK), intent(in)                 :: d
    integer(IK), intent(inout)              :: X(d, *)
    integer(IK)                             :: i, j
!
    if (.not. ALLOCATED(this%p)) return
!
    do concurrent(i=1:SIZE(this%p, 2))
      do concurrent(j=1:this%p(3, i))
        block
          integer(IK) :: b
          b = this%p(1, i) + this%p(2, i) * (j - 1)
          call cyclic_swap_int(d, this%p(2, i), this%q(b), X)
        end block
      end do
    end do
!
  end subroutine group_permutation_swap_int_l2
!
  pure subroutine group_permutation_swap_real_l1(this, X)
    class(group_permutation), intent(inout) :: this
    real(RK), intent(inout)                 :: X(*)
    integer(IK)                             :: i, j
!
    if (.not. ALLOCATED(this%p)) return
!
    do concurrent(i=1:SIZE(this%p, 2))
      do concurrent(j=1:this%p(3, i))
        block
          integer(IK) :: b
          b = this%p(1, i) + this%p(2, i) * (j - 1)
          call cyclic_swap_real(1, this%p(2, i), this%q(b), X)
        end block
      end do
    end do
!
  end subroutine group_permutation_swap_real_l1
!
  pure subroutine group_permutation_swap_real_l2(this, d, X)
    class(group_permutation), intent(inout) :: this
    integer(IK), intent(in)                 :: d
    real(RK), intent(inout)                 :: X(d, *)
    integer(IK)                             :: i, j
!
    if (.not. ALLOCATED(this%p)) return
!
    do concurrent(i=1:SIZE(this%p, 2))
      do concurrent(j=1:this%p(3, i))
        block
          integer(IK) :: b
          b = this%p(1, i) + this%p(2, i) * (j - 1)
          call cyclic_swap_real(d, this%p(2, i), this%q(b), X)
        end block
      end do
    end do
!
  end subroutine group_permutation_swap_real_l2
!
  pure elemental function group_permutation_nfree(this) result(res)
    class(group_permutation), intent(in) :: this
    integer(IK)                          :: res
    if (ALLOCATED(this%q)) then
      res = SIZE(this%q)
    else
      res = 0
    end if
  end function group_permutation_nfree
!
  pure function group_permutation_free_indices(this) result(res)
    class(group_permutation), intent(in) :: this
    integer(IK)                          :: res(this%nfree())
    integer(IK)                          :: i, n
    if (ALLOCATED(this%q)) then
      n = SIZE(res)
      do concurrent(i=1:n)
        res(i) = this%q(i)
      end do
      call qsi(n, res)
    end if
  end function group_permutation_free_indices
!
  pure elemental subroutine group_permutation_clear(this)
    class(group_permutation), intent(inout) :: this
    if (ALLOCATED(this%p)) deallocate (this%p)
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine group_permutation_clear
!
  pure elemental subroutine group_permutation_destroy(this)
    type(group_permutation), intent(inout) :: this
    call this%clear()
  end subroutine group_permutation_destroy
!
!!!
!
  pure subroutine cyclic_swap_int(d, s, q, X)
    integer(IK), intent(in)    :: d, s, q(*)
    integer(IK), intent(inout) :: X(d, *)
    integer(IK)                :: T(d)
    integer(IK)                :: i, j
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
  end subroutine cyclic_swap_int
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
  pure subroutine uniq(n, a, u, res)
    integer(IK), intent(in)    :: n, a(*)
    integer(IK), intent(inout) :: u, res(*)
    integer(IK)                :: i
!
    u = 0
    if (n < 1) return
!
    do concurrent(i=1:n)
      res(i) = a(i)
    enddo
    call qsi(n, res)
!
    u = 1
    res(n + u) = 1
    do i = 1, n - 1
      if (res(i) == res(i + 1)) then
        res(n + u) = res(n + u) + 1
        cycle
      endif
      if (res(i) > 1) u = u + 1
      res(u) = res(i + 1)
      res(n + u) = 1
    end do
!
    do i = 1, u
      res(u + i) = res(n + i)
    end do
!
  end subroutine uniq
!
  pure recursive subroutine qsi(n, s)
    integer(IK), intent(in)    :: n
    integer(IK), intent(inout) :: s(n)
    integer(IK)                :: i, j, p, t
!
    j = SIZE(s); i = 1; p = j / 2
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
    if (j + 1 < SIZE(s)) call qsi(n - j, s(j + 1:))
!
  end subroutine qsi
!
end module mod_group_permutation
