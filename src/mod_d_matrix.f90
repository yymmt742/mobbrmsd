module mod_d_matrix
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_molecular_rotation
  use mod_mol_block
  use mod_Kabsch
  use mod_Hungarian
  implicit none
  private
  public :: d_matrix
  public :: d_matrix_memsize
  public :: d_matrix_eval
  public :: d_matrix_partial_lower_bound
!
  type d_matrix
    private
    sequence
    integer(IK) :: d, s, m, n
    integer(IK) :: dd, dm, dds, ddsn
    integer(IK) :: z, c, h
  end type d_matrix
!
  interface d_matrix
    module procedure d_matrix_new
  end interface d_matrix
!
  interface
    include 'dgemm.h'
  end interface
!
contains
!
!| generator
  pure elemental function d_matrix_new(p, d, s, b) result(res)
    integer(IK), intent(in)     :: p, d, s
    type(mol_block), intent(in) :: b
    type(d_matrix)              :: res
    res%d = MAX(d, 1)
    res%s = MAX(s, 1)
    res%m = b%m
    res%n = b%g
    res%dd = res%d * res%d
    res%dm = res%d * res%m
    res%dds = res%dd * res%s
    res%ddsn = res%dds * res%n
    res%h  = p
    res%z  = res%h + 1
    res%c  = res%z + res%n * res%n
  end function d_matrix_new
!
  pure elemental function d_matrix_memsize(a) result(res)
    type(d_matrix), intent(in) :: a
    integer(IK)                :: res
    res = (a%dds + a%s + 1) * a%n * a%n + 1
  end function d_matrix_memsize
!
  pure subroutine d_matrix_eval(a, rot, X, Y, W)
    type(d_matrix), intent(in)            :: a
    class(molecular_rotation), intent(in) :: rot
    real(RK), intent(in)                  :: X(*)
    real(RK), intent(in)                  :: Y(*)
    real(RK), intent(inout)               :: W(*)
!
    call eval(a%d, a%s, a%m, a%n, a%dd, a%dm, rot, X, Y, W(a%h), W(a%z), W(a%c))
!
  end subroutine d_matrix_eval
!
  pure subroutine eval(d, s, m, n, dd, dm, r, X, Y, H, Z, C)
    integer(IK), intent(in)              :: d, s, m, n
    integer(IK), intent(in)              :: dd, dm
    type(molecular_rotation), intent(in) :: r
    real(RK), intent(in)                 :: X(d, m, n)
    real(RK), intent(in)                 :: Y(d, m, n)
    real(RK), intent(inout)              :: H(1)
    real(RK), intent(inout)              :: Z(n, n)
    real(RK), intent(inout)              :: C(dd + 1, s, n, n)
    integer(IK), parameter               :: ib = 1
    integer(IK), parameter               :: it = 2
    integer(IK), parameter               :: ih = 3
    integer(IK), parameter               :: ix = 4
    integer(IK)                          :: iy, ic, ir, iw
    integer(IK)                          :: j, k, nw
!
    if (n < 1) return
!
    nw = nwork(d, m)
    iy = ix + dm
    ic = iy + dm
    ir = ic + dd
    iw = ir + dd
!
    do concurrent(j=1:n, k=1:n)
      block
        integer(IK) :: i
        real(RK)    :: W(nw)
!
        call copy(dm, X(1, 1, j), W(ix))
        call copy(dm, Y(1, 1, k), W(iy))
!
        call calc_lb(d, m, dd, dm, w(ix), W(ih), W(it), W(ic), W(ir), W(iw), C(1, 1, j, k))
        w(ib) = w(it)
!
        do i = 2, s
          call r%swap(d, W(iy), i - 1)
          call calc_lb(d, m, dd, dm, W(ix), W(ih), W(it), W(ic), W(ir), W(iw), C(1, i, j, k))
          w(ib) = MIN(w(ib), w(it))
          call r%reverse(d, W(iy), i - 1)
        end do
!
        Z(j, k) = w(ib)
!
      end block
    enddo
!
!!! inquire work memory size
    call Hungarian(-n, Z, H)
    nw = NINT(H(1), IK)
!
    block
      real(RK)    :: W(nw)
      call Hungarian(n, Z, W)
      H(1) = W(1)
    end block
!
  contains
!
    pure elemental function nwork(d, dm) result(res)
      integer(IK), intent(in) :: d, dm
      integer(IK)             :: res
      res = 3 + 2 * d * d + dm + dm + Kabsch_worksize(d)
    end function nwork
!
  end subroutine eval
!
  pure subroutine calc_lb(d, m, dd, dm, XY, H, T, C, R, W, CT)
    integer(IK), intent(in) :: d, m, dd, dm
    real(RK), intent(in)    :: XY(dm, *)
    real(RK), intent(inout) :: H, T, C(d, d), R(d, d), W(*), CT(*)
!!! trace of self correlation matrix
    H = ddot(dm + dm, XY, XY) ! tr[X, X^t] + tr[Y, Y^t]
!!! get correlation matrix C = Y^t@X and optimal rotation R^t
    call DGEMM('N', 'T', d, d, m, ONE, XY(1, 2), d, XY(1, 1), d, ZERO, C, d)
    call Kabsch(d, C, R, W)
!!! get squared displacement
    T = ddot(dd, C, R) ! tr[C, R^t]
    T = T + T
    T = H - T          ! SD = tr[X^tX] + tr[Y^tY] - 2*tr[C, R]
!!! summarize to memory
    CT(1) = H
    call copy(dd, C, CT(2))
  end subroutine calc_lb
!
  pure function d_matrix_partial_lower_bound(a, np, list, W) result(res)
    type(d_matrix), intent(in) :: a
    integer(IK), intent(in)    :: np, list(*)
    real(RK), intent(in)       :: W(*)
    real(RK)                   :: res
    integer(IK)                :: nw, pp, iw
!
    if (np < 1 .or. a%n == np) then
      res = W(a%h)
    elseif (a%n < np) then
      res = ZERO
    else
      nw = nwork(np)
      pp = np * np
      iw = 1 + pp
      block
        real(RK) :: T(nw)
        call pack_Z(a%n, np, list, W(a%z), T(1))
        call Hungarian(np, T(1), T(iw))
        res = T(iw)
      end block
    end if
!
  contains
!
    pure elemental function nwork(np) result(res)
      integer(IK), intent(in) :: np
      real(RK)                :: W(1)
      integer(IK)             :: res
      call Hungarian(-np, W, W)
      res = np * np + NINT(W(1), IK)
    end function nwork
!
  end function d_matrix_partial_lower_bound
!
  pure subroutine pack_Z(n, np, list, Z, T)
    integer(IK), intent(in) :: n, np, list(np)
    real(RK), intent(in)    :: Z(n, n)
    real(RK), intent(inout) :: T(np, np)
    integer(IK)             :: i, j
!
    do concurrent(i=1:np, j=1:np)
      T(i, j) = Z(list(i), n - np + j)
    end do
!
  end subroutine pack_Z
!
  pure function ddot(d, X, Y) result(res)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: X(*), Y(*)
    real(RK)                :: res
    if(d==1)then
      res = X(1) * Y(1)
    elseif(d==4)then
      res = X(1) * Y(1) + X(2) * Y(2) + &
        &   X(3) * Y(3) + X(4) * Y(4)
    elseif(d==9)then
      res = X(1) * Y(1) + X(2) * Y(2) + X(3) * Y(3) +&
        &   X(4) * Y(4) + X(5) * Y(5) + X(6) * Y(6) +&
        &   X(7) * Y(7) + X(8) * Y(8) + X(9) * Y(9)
    else
      block
        integer(IK) :: i
        res = ZERO
        do i = 1, d
          res = res + X(i) * Y(i)
        end do
      end block
    endif
  end function ddot
!
  pure subroutine copy(d, source, dest)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: source(*)
    real(RK), intent(inout) :: dest(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dest(i) = source(i)
    end do
  end subroutine copy
!
end module mod_d_matrix
