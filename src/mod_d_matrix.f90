module mod_d_matrix
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_molecular_rotation
  use mod_mol_block
  use mod_Kabsch
  use mod_Hungarian
  implicit none
  private
  public :: d_matrix
  public :: d_matrix_list
  public :: d_matrix_memsize
  public :: d_matrix_eval
  public :: d_matrix_partial_eval
!
  type d_matrix
    private
    sequence
    integer(IK) :: d, s, m, n
    integer(IK) :: dd, nn, dm, cb
    integer(IK) :: x, z, c, nk
  end type d_matrix
!
  interface d_matrix
    module procedure d_matrix_new
  end interface d_matrix
!
  type d_matrix_list
    integer(IK)                           :: d = 0
    integer(IK)                           :: n = 0
    type(d_matrix), allocatable           :: m(:)
    type(molecular_rotation), allocatable :: r(:)
  contains
    procedure :: eval => d_matrix_list_eval
    final     :: d_matrix_list_destroy
  end type d_matrix_list
!
  interface d_matrix_list
    module procedure d_matrix_list_new
  end interface d_matrix_list
!
  interface
    include 'dgemm.h'
  end interface
!
contains
!
!| generator
  pure elemental function d_matrix_new(p, q, d, s, b) result(res)
    integer(IK), intent(in)     :: p, q, d, s
    type(mol_block), intent(in) :: b
    type(d_matrix)              :: res
    res%d = MAX(d, 1)
    res%s = MAX(s, 1)
    res%m = b%m
    res%n = b%g
    res%dd = res%d * res%d
    res%nn = res%n * res%n
    res%dm = res%d * res%m
    res%cb = res%dd * res%s + 1
    res%x  = q
    res%z  = p
    res%c  = res%z + res%nn
    res%nk = Kabsch_worksize(d)
  end function d_matrix_new
!
  pure elemental function d_matrix_memsize(a) result(res)
    type(d_matrix), intent(in) :: a
    integer(IK)                :: res
    res = (a%cb + 1) * a%nn
  end function d_matrix_memsize
!
  pure subroutine d_matrix_eval(a, rot, X, Y, W)
    type(d_matrix), intent(in)            :: a
    class(molecular_rotation), intent(in) :: rot
    real(RK), intent(in)                  :: X(*)
    real(RK), intent(in)                  :: Y(*)
    real(RK), intent(inout)               :: W(*)
!
    call eval(a%d, a%s, a%m, a%n, a%dd, a%dm, a%cb, a%nk, rot, X(a%x), Y(a%x), W(a%z), W(a%c))
!
  end subroutine d_matrix_eval
!
  pure subroutine eval(d, s, m, n, dd, dm, cb, nk, r, X, Y, Z, C)
    integer(IK), intent(in)              :: d, s, m, n
    integer(IK), intent(in)              :: dd, dm, cb, nk
    type(molecular_rotation), intent(in) :: r
    real(RK), intent(in)                 :: X(d, m, n)
    real(RK), intent(in)                 :: Y(d, m, n)
    real(RK), intent(inout)              :: Z(n, n)
    real(RK), intent(inout)              :: C(cb, n, n)
    integer(IK), parameter               :: ib = 1
    integer(IK), parameter               :: it = 2
    integer(IK), parameter               :: ih = 3
    integer(IK), parameter               :: ix = 4
    integer(IK)                          :: iy, ic, ir, iw
    integer(IK)                          :: j, k, nw
!
    if (n < 1) return
!
    nw = 3 + 2 * d * d + dm + dm + nk
    iy = ix + dm
    ic = iy + dm
    ir = ic + dd
    iw = ir + dd
!
    do concurrent(j=1:n, k=1:n)
      block
        integer(IK) :: i, ip
        real(RK)    :: W(nw)
!
        call copy(dm, X(1, 1, j), W(ix))
        call copy(dm, Y(1, 1, k), W(iy))
!
!!!     trace of self correlation matrix
        w(ih) = ddot(dm + dm, W(ix), W(ix)) ! tr[X, X^t] + tr[Y, Y^t]
        C(1, j, k) = w(ih)
!
        ip = 2
        call calc_lb(d, m, dd, dm, w(ih), w(ix), W(it), W(ic), W(ir), W(iw))
        call copy(dd, W(ic), C(ip, j, k))
        w(ib) = w(it)
!
        do i = 1, s - 1
          call r%swap(d, W(iy), i)
          ip = ip + dd
          call calc_lb(d, m, dd, dm, W(ih), W(ix), W(it), W(ic), W(ir), W(iw))
          call copy(dd, W(ic), C(ip, j, k))
          w(ib) = MIN(w(ib), w(it))
          call r%reverse(d, W(iy), i)
        end do
!
        Z(j, k) = w(ib)
!
      end block
    enddo
!
  contains
!
    pure subroutine calc_lb(d, m, dd, dm, H, XY, T, C, R, W)
      integer(IK), intent(in) :: d, m, dd, dm
      real(RK), intent(in)    :: H, XY(dm, *)
      real(RK), intent(inout) :: T, C(d, d), R(d, d), W(*)
!!!   get correlation matrix C = Y^t@X and optimal rotation R^t
      call DGEMM('N', 'T', d, d, m, ONE, XY(1, 2), d, XY(1, 1), d, ZERO, C, d)
      call Kabsch(d, C, R, W)
!!!   get squared displacement
      T = ddot(dd, C, R) ! tr[C, R^t]
      T = T + T
      T = H - T
    end subroutine calc_lb
!
  end subroutine eval
!
  pure subroutine d_matrix_partial_eval(a, p, iprm, isym, ires, W, LF, LB, H, C, R)
    type(d_matrix), intent(in)        :: a
    integer(IK), intent(in)           :: p, iprm, isym, ires(*)
    real(RK), intent(in)              :: W(*)
    real(RK), intent(inout)           :: LF, LB, H, C(*)
    real(RK), intent(inout), optional :: R(*)
    integer(IK)                       :: ih, ic, nw
!
    if (p < 0 .or. a%n < p) return
    if (iprm < 1 .or. a%n < iprm) return
    if (isym < 1 .or. a%s < isym) return
!
    ih = a%c + (iprm - 1) * a%cb + (p - 1) * a%cb * a%n
    ic = ih + 1 + a%dd * (isym - 1)
    nw = 1 + a%dd * 2 + a%nk
!
    if (p == 0) then
      LF = ZERO
      if (PRESENT(R)) call eye(a%d, R)
    else
      call partial_eval(a%d, a%s, a%n, a%dd, nw, p, iprm, isym, W(ih), W(ic), LF, H, C, R)
    end if
!
    if (p == a%n) then
      LB = ZERO
    else
      call residue_eval(a, p, ires, W, LB)
    endif
!
  end subroutine d_matrix_partial_eval
!
  pure subroutine partial_eval(d, s, n, dd, nw, p, iprm, isym, H, C, LF, HP, CP, R)
    integer(IK), intent(in) :: d, s, n, dd, nw
    integer(IK), intent(in) :: p, iprm, isym
    real(RK), intent(in)    :: H, C(*)
    real(RK), intent(inout) :: LF, HP, CP(*)
    real(RK), intent(inout), optional :: R(*)
    real(RK)                :: W(nw)
    integer(IK), parameter  :: it = 1
    integer(IK), parameter  :: ic = 2
    integer(IK)             :: ir, iw
!
    ir = ic + dd
    iw = ir + dd
!
!!! update H and C
    w(it) = HP + H
    HP = w(it)
    call add(dd, CP, C, W(ic))
    call copy(dd, w(ic), CP)
!
!!! get correlation matrix C = Y^t@X and optimal rotation R^t
    call Kabsch(d, w(ic), w(ir), W(iw))
!!! get squared displacement
    W(it) = ddot(dd, W(ic), W(ir)) ! tr[C, R^t]
    W(it) = W(it) + W(it)
    LF = HP - w(it)
!
!!! summarize to memory
    if (PRESENT(R)) call copy(dd, w(ir), R)
!
  end subroutine partial_eval
!
  pure subroutine residue_eval(a, p, ires, W, res)
    type(d_matrix), intent(in) :: a
    integer(IK), intent(in)    :: p, ires(*)
    real(RK), intent(in)       :: W(*)
    real(RK), intent(inout)    :: res
    integer(IK)                :: nw, nr, rr
!
    if (p < 0 .or. a%n < p) then
      res = ZERO
    elseif (p == 0) then
      nw = nwork(a%n)
      block
        real(RK) :: T(nw)
        call Hungarian(a%n, W(a%z), T(1))
        res = T(1)
      end block
    else
      nr = a%n - p
      rr = nr * nr
      nw = rr + nwork(nr)
      block
        integer(IK) :: iw, iz
        real(RK)    :: T(nw)
        iw = 1 + rr
        iz = a%z + a%n * (p - 1)
        call pack_Z(a%n, nr, ires, W(a%z), T(1))
        call Hungarian(nr, T(1), T(iw))
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
      res = NINT(W(1), IK)
    end function nwork
!
    pure subroutine pack_Z(n, nr, ires, Z, T)
      integer(IK), intent(in) :: n, nr, ires(nr)
      real(RK), intent(in)    :: Z(n, *)
      real(RK), intent(inout) :: T(nr, nr)
      integer(IK)             :: i, j
      do concurrent(i=1:nr, j=1:nr)
        T(i, j) = Z(ires(i), n - nr + j)
      end do
    end subroutine pack_Z
!
  end subroutine residue_eval
!
!!! d_matrix_list
!
  pure function d_matrix_list_new(n, d, b, r) result(res)
    integer(IK), intent(in)              :: n, d
    type(mol_block), intent(in)          :: b(n)
    type(molecular_rotation), intent(in) :: r(n)
    type(d_matrix_list)                  :: res
    integer(IK)                          :: i, p, q, s
!
    res%n = MAX(n, 0)
    res%d = MAX(d, 1)
!
    allocate(res%m(res%n))
    allocate(res%r(res%n))
!
    do concurrent(i=1:res%n)
      res%r(i) = r(i)
    end do
!
    p = 1
    q = 1
!
    do concurrent(i=1:res%n)
      s = r(i)%n_sym()
      res%m(i) = d_matrix(p, q, d, s, b(i))
      p = p + d_matrix_memsize(res%m(i))
      q = q + d * b(i)%m * b(i)%n
    end do
!
  end function d_matrix_list_new
!
  pure subroutine d_matrix_list_eval(this, X, Y, W)
    class(d_matrix_list), intent(in) :: this
    real(RK), intent(in)             :: X(*), Y(*)
    real(RK), intent(inout)          :: W(*)
    integer(IK)                      :: i
    do concurrent(i=1:this%n)
      call d_matrix_eval(this%m(i), this%r(i), X, Y, W)
    enddo
  end subroutine d_matrix_list_eval
!
  pure elemental subroutine d_matrix_list_destroy(this)
    type(d_matrix_list), intent(inout) :: this
    this%n = 0
    this%d = 0
    if (ALLOCATED(this%m)) deallocate (this%m)
    if (ALLOCATED(this%r)) deallocate (this%r)
  end subroutine d_matrix_list_destroy
!
!!! util
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
  pure subroutine add(d, A, B, C)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: A(*), B(*)
    real(RK), intent(inout) :: C(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      C(i) = A(i) + B(i)
    end do
  end subroutine add
!
  pure subroutine eye(d, x)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: x(d, *)
    integer(IK)             :: i, j
    do concurrent(j=1:d, i=1:d)
      x(i, j) = MERGE(ONE, ZERO, i == j)
    end do
  end subroutine eye
!
end module mod_d_matrix
