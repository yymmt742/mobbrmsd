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
    integer(IK), public :: s, m, n, g
    integer(IK)         :: d, dd, gg, dm, cb
    integer(IK)         :: x, z, c, nk
  end type d_matrix
!
  interface d_matrix
    module procedure d_matrix_new
  end interface d_matrix
!
  type d_matrix_list
    integer(IK)                           :: h = 0
    integer(IK)                           :: v = 0
    integer(IK)                           :: c = 0
    integer(IK)                           :: o = 0
    integer(IK)                           :: d = 0
    integer(IK)                           :: l = 0
    integer(IK)                           :: dd = 0
    type(d_matrix), allocatable           :: m(:)
  contains
    procedure :: n_depth      => d_matrix_list_n_depth
    procedure :: eval         => d_matrix_list_eval
    procedure :: partial_eval => d_matrix_list_partial_eval
    procedure :: memsize      => d_matrix_list_memsize
    procedure :: clear        => d_matrix_list_clear
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
  pure elemental function d_matrix_new(p, d, b) result(res)
    integer(IK), intent(in)     :: p, d
    type(mol_block), intent(in) :: b
    type(d_matrix)              :: res
    res%d = MAX(d, 1)
    res%s = MAX(b%s, 1)
    res%m = MAX(b%m, 0)
    res%n = MAX(b%n, 0)
    if(b%g==1.and.res%s==1)then
      res%g = 0
    else
      res%g = MIN(res%n, MAX(b%g, 0))
    endif
    res%dd = res%d * res%d
    res%gg = res%g * res%g
    res%dm = res%d * res%m
    res%cb = res%dd * res%s + 1
    res%x  = b%p
    res%z  = p
    res%c  = res%z + res%gg
    res%nk = Kabsch_worksize(d)
  end function d_matrix_new
!
  pure elemental function d_matrix_memsize(a) result(res)
    type(d_matrix), intent(in) :: a
    integer(IK)                :: res
    res = (a%cb + 1) * a%gg
  end function d_matrix_memsize
!
  pure subroutine d_matrix_eval(a, rot, X, Y, W)
    type(d_matrix), intent(in)            :: a
    class(molecular_rotation), intent(in) :: rot
    real(RK), intent(in)                  :: X(*)
    real(RK), intent(in)                  :: Y(*)
    real(RK), intent(inout)               :: W(*)
!
    call eval(a%d, a%s, a%m, a%g, a%dd, a%dm, a%cb, a%nk, rot, X(a%x), Y(a%x), W(a%z), W(a%c))
!
  end subroutine d_matrix_eval
!
  pure subroutine eval(d, s, m, g, dd, dm, cb, nk, r, X, Y, Z, C)
    integer(IK), intent(in)              :: d, s, m, g
    integer(IK), intent(in)              :: dd, dm, cb, nk
    type(molecular_rotation), intent(in) :: r
    real(RK), intent(in)                 :: X(d, m, g)
    real(RK), intent(in)                 :: Y(d, m, g)
    real(RK), intent(inout)              :: Z(g, g)
    real(RK), intent(inout)              :: C(cb, g, g)
    integer(IK), parameter               :: ib = 1
    integer(IK), parameter               :: it = 2
    integer(IK), parameter               :: ih = 3
    integer(IK), parameter               :: ix = 4
    integer(IK)                          :: iy, ic, ir, iw
    integer(IK)                          :: j, k, nw
!
    if (g < 1) return
!
    nw = 3 + 2 * d * d + dm + dm + nk
    iy = ix + dm
    ic = iy + dm
    ir = ic + dd
    iw = ir + dd
!
    do concurrent(j=1:g, k=1:g)
      block
        integer(IK) :: i, ip
        real(RK)    :: W(nw)
!
        call copy(dm, X(1, 1, j), W(ix))
        call copy(dm, Y(1, 1, k), W(iy))
!
!!!     trace of self correlation matrix
        call ddot(dm + dm, W(ix), W(ix), w(ih)) ! tr[X, X^t] + tr[Y, Y^t]
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
      call ddot(dd, C, R, T) ! tr[C, R^t]
      T = T + T
      T = H - T
    end subroutine calc_lb
!
  end subroutine eval
!
  pure subroutine d_matrix_partial_eval(a, p, iprm, isym, ires, W, LT, H, C, LF, LB, R)
    type(d_matrix), intent(in)        :: a
    integer(IK), intent(in)           :: p, iprm, isym, ires(*)
    real(RK), intent(in)              :: W(*)
    real(RK), intent(inout)           :: LT, H, C(*)
    real(RK), intent(inout), optional :: LF, LB, R(*)
!
    LT = ZERO
    if (PRESENT(LF)) LF = ZERO
    if (PRESENT(LB)) LB = ZERO
    if (PRESENT(R)) call eye(a%d, R)
!
    if (p < 0 .or. a%g < p) return
    if (p < a%g) call setminus_eval(a, p, ires, W, LT)
    if (PRESENT(LB)) LB = LT
!
    if (p == 0) return
    if (iprm < 1 .or. a%g < iprm) return
    if (isym < 0 .or. a%s <= isym) return
!
    block
      integer(IK) :: ih, ic, nw
      ih = a%c + (iprm - 1) * a%cb + (p - 1) * a%cb * a%g
      ic = ih + 1 + a%dd * isym
      nw = 1 + a%dd * 2 + a%nk
      if (PRESENT(LF)) then
        call partial_eval(a%d, a%s, a%g, a%dd, nw, p, W(ih), W(ic), LF, H, C, R)
        LT = LT + LF
      else
        call partial_eval(a%d, a%s, a%g, a%dd, nw, p, W(ih), W(ic), LT, H, C, R)
      end if
    end block
!
  contains
!
    pure subroutine setminus_eval(a, p, ires, W, res)
      type(d_matrix), intent(in) :: a
      integer(IK), intent(in)    :: p, ires(*)
      real(RK), intent(in)       :: W(*)
      real(RK), intent(inout)    :: res
      integer(IK)                :: nw, nr, rr
!
      if (p < 0 .or. a%g < p) then
        res = ZERO
      elseif (p == 0) then
        nw = nwork(a%g)
        block
          real(RK) :: T(nw)
          call Hungarian(a%g, W(a%z), T(1))
          res = T(1)
        end block
      else
        nr = a%g - p
        rr = nr * nr
        nw = rr + nwork(nr)
        block
          integer(IK) :: iw, iz
          real(RK)    :: T(nw)
          iw = 1 + rr
          iz = a%z + a%g * (p - 1)
          call pack_Z(a%g, nr, ires, W(a%z), T(1))
          call Hungarian(nr, T(1), T(iw))
          res = T(iw)
        end block
      end if
!
    end subroutine setminus_eval
!
    pure subroutine partial_eval(d, s, n, dd, nw, p, H, C, LF, HP, CP, R)
      integer(IK), intent(in) :: d, s, n, dd, nw
      integer(IK), intent(in) :: p
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
      call ddot(dd, W(ic), W(ir), W(it)) ! tr[C, R^t]
      W(it) = W(it) + W(it)
      LF = LF + HP - w(it)
!
!!! summarize to memory
      if (PRESENT(R)) call copy(dd, w(ir), R)
!
    end subroutine partial_eval
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
  end subroutine d_matrix_partial_eval
!
!!! d_matrix_list
!
  pure function d_matrix_list_new(b, p) result(res)
    type(mol_block_list), intent(in)     :: b
    integer(IK), intent(in)              :: p
    type(d_matrix_list)                  :: res
    integer(IK)                          :: i, ip
!
    res%d = b%nspatial()
    res%l = b%nspecies()
    res%dd = res%d**2
    res%h = p              ! H(1)
    res%v = res%h + 1      ! V(1)
    res%c = res%v + 1      ! C(d*d)
    res%o = res%c + res%dd ! O(L+1)
    ip = res%o + res%l + 1
!
    allocate (res%m(res%l))
!
    do i = 1, res%l
      res%m(i) = d_matrix(ip, res%d, b%b(i))
      ip = ip + d_matrix_memsize(res%m(i))
    end do
!
  end function d_matrix_list_new
!
  pure function d_matrix_list_memsize(this) result(res)
    class(d_matrix_list), intent(in) :: this
    integer(IK)                      :: res
    if (ALLOCATED(this%m)) then
      res = SUM(d_matrix_memsize(this%m)) + (this%l + 1) + this%d**2 + 2
    else
      res = 0
    end if
  end function d_matrix_list_memsize
!
  pure elemental function d_matrix_list_n_depth(this) result(res)
    class(d_matrix_list), intent(in) :: this
    integer(IK)                      :: res
    if (ALLOCATED(this%m)) then
      res = SUM(this%m%g)
    else
      res = 0
    end if
  end function d_matrix_list_n_depth
!
  pure subroutine d_matrix_list_eval(this, rot, X, Y, W)
    class(d_matrix_list), intent(in)     :: this
    type(molecular_rotation), intent(in) :: rot(*)
    real(RK), intent(in)                 :: X(*), Y(*)
    real(RK), intent(inout)              :: W(*)
    integer(IK)                          :: i
!
    if (.not. ALLOCATED(this%m)) return
!
    call fixpoints_eval(this%d, this%l, this%m, X, Y, W(this%h), W(this%v), W(this%c))
!
    do concurrent(i=1:this%l)
      call d_matrix_eval(this%m(i), rot(i), X, Y, W)
      block
        integer(IK) :: j, ires(this%m(i)%g)
        do concurrent(j=1:this%m(i)%g)
          ires(j) = j
        end do
        j = this%o + i - 1
        call d_matrix_partial_eval(this%m(i), 0, 0, 0, ires, W, W(j), W(j), W(j))
      end block
    end do
!
    do concurrent(i=this%l - 1:1:-1)
      W(this%o + i - 1) = W(this%o + i - 1) + W(this%o + i)
    end do
!
    W(this%o + this%l) = ZERO
!
  contains
!
    pure subroutine fixpoints_eval(d, l, m, X, Y, H, V, C)
      integer(IK), intent(in)    :: d, l
      type(d_matrix), intent(in) :: m(l)
      real(RK), intent(in)       :: X(*), Y(*)
      real(RK), intent(inout)    :: H, V, C(*)
      integer(IK)                :: t(l), p(l), q(l)
      integer(IK), parameter     :: ih = 1
      integer(IK), parameter     :: iv = 2
      integer(IK), parameter     :: ic = 3
      integer(IK)                :: i, dd, nk, ir, iw, ix, iy, mn, dmn, nw
!
      do concurrent(i=1:l)
        t(i) = m(i)%m * (m(i)%n - m(i)%g)
      end do
!
      nk = Kabsch_worksize(d)
      dd = d * d
      mn = SUM(t)
      dmn = d * mn
!
      ir = ic + dd
      iw = ir + dd
      ix = ic + dd
      iy = ix + dmn
      nw = 2 + dd + MAX(dmn + dmn, dd + nk)
!
      do concurrent(i=1:l)
        t(i) = t(i) * d
      end do
      p(1) = 0
      do i = 2, l
        p(i) = p(i - 1) + t(i - 1)
      end do
      do concurrent(i=1:l)
        q(i) = m(i)%x + m(i)%dm * m(i)%g
      end do
!
      block
        real(RK) :: W(nw)
!
        do concurrent(i=1:l)
          block
            integer(IK) :: px
            px = p(i) + ix
            call copy(t(i), X(q(i)), W(px))
          end block
        end do
        do concurrent(i=1:l)
          block
            integer(IK) :: py
            py = p(i) + iy
            call copy(t(i), Y(q(i)), W(py))
          end block
        end do
!
        call ddot(dmn + dmn, W(ix), W(ix), W(ih))
!
        if (mn > 0) then
          call DGEMM('N', 'T', d, d, mn, ONE, W(iy), d, W(ix), d, ZERO, W(ic), d)
          call Kabsch(d, w(ic), w(ir), W(iw))
          call ddot(dd, w(ic), w(ir), w(iv))
        else
          call zfill(dd, W(ic))
          w(iv) = ZERO
        end if
!
        H = W(ih)
        V = W(ih) - W(iv) - W(iv)
        call copy(dd, W(ic), C)
!
      end block
!
    end subroutine fixpoints_eval
!
  end subroutine d_matrix_list_eval
!
  pure subroutine d_matrix_list_partial_eval(this, p, perm, iprm, isym, W, LT, H, C, LF, LB, R)
    class(d_matrix_list), intent(in)  :: this
    integer(IK), intent(in)           :: p, perm(*), iprm, isym
    real(RK), intent(in)              :: W(*)
    real(RK), intent(inout)           :: LT, H, C(*)
    real(RK), intent(inout), optional :: LF, LB, R(*)
    integer(IK)                       :: ispc, iofs, nres
!
    call p_index(this, p, ispc, iofs)
    if (ispc < 1) return
!
    nres = this%m(ispc)%g - iofs
    block
      integer(IK) :: i, jper, ires(nres)
      jper = perm(p + iprm - 1)
      do concurrent(i=1:iprm - 1)
        ires(i) = perm(p + i - 1)
      end do
      do concurrent(i=iprm:nres)
        ires(i) = perm(p + i)
      end do
      call d_matrix_partial_eval(this%m(ispc), iofs, jper, isym, ires, W, LT, H, C, LF, LB, R)
    end block
!
  contains
!
    pure elemental subroutine p_index(this, p, ispc, iofs)
      type(d_matrix_list), intent(in) :: this
      integer(IK), intent(in)         :: p
      integer(IK), intent(inout)      :: ispc, iofs
      integer(IK)                     :: i, q, r
!
      ispc = 0
      iofs = 0
      r = p
      q = 1
!
      do i = 1, this%l
        if (q > p) return
        ispc = i
        iofs = r
        r = r - this%m(i)%g
        q = q + this%m(i)%g
      end do
!
      if (q <= p) ispc = 0
!
    end subroutine p_index
!
  end subroutine d_matrix_list_partial_eval
!
  pure elemental subroutine d_matrix_list_clear(this)
    class(d_matrix_list), intent(inout) :: this
    this%l = 0
    this%d = 0
    if (ALLOCATED(this%m)) deallocate (this%m)
  end subroutine d_matrix_list_clear
!
  pure elemental subroutine d_matrix_list_destroy(this)
    type(d_matrix_list), intent(inout) :: this
    call d_matrix_list_clear(this)
  end subroutine d_matrix_list_destroy
!
!!! util
!
  pure subroutine ddot(d, X, Y, res)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: X(*), Y(*)
    real(RK), intent(inout) :: res
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
  end subroutine ddot
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
  pure subroutine zfill(d, x)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: x(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      x(i) = ZERO
    end do
  end subroutine zfill
!
end module mod_d_matrix
