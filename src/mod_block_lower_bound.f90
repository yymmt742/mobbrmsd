!| Calculate min_{R,P,Q1,...,Q_n} |X-R@Y@f(P)@g(Q1,Q2,...,Qn)|^2 for R\in\R^{d,d}, P\in \R^{n,n} and Qi\in \R^{m,m}.
!  f(P) = I x P ( I is the identity and \in R^{m,m} ).
!  g(Q1,Q2,...,Qn) = Diag(Q1 Q2 ... Qn).
!  Here, R@R^T=I, det(R)=1, P@P^T=I, and Q@Q^T=I are satisfied.
module mod_block_lower_bound
  use mod_params, only: IK, RK, ONE => RONE, FOUR => RFOUR, ZERO => RZERO, RHUGE
  use mod_optarg
  use mod_rmsd
  use mod_det
  use mod_cov
  use mod_Kabsch
  use mod_Procrustes
  implicit none
  private
  public :: block_lower_bound_worksize, block_lower_bound
!
  type mol_block
    sequence
    integer(IK) :: m = 1
    integer(IK) :: n = 1
    integer(IK) :: g = 1
    integer(IK) :: f = 1
  end type mol_block
!
  type mol_block_
    private
    sequence
    integer(IK) :: m
    integer(IK) :: n
    integer(IK) :: g
    integer(IK) :: f
    integer(IK) :: mg
    integer(IK) :: mn
    integer(IK) :: df
    integer(IK) :: dm
    integer(IK) :: ff
    integer(IK) :: gg
    integer(IK) :: dmn
    integer(IK) :: dmg
    integer(IK) :: m_f
    integer(IK) :: ffg
    integer(IK) :: ffgg
    integer(IK) :: dxm_dxf
    integer(IK) :: proc_f
    integer(IK) :: blksiz
    integer(IK) :: ix
    integer(IK) :: iy
    integer(IK) :: ip
    integer(IK) :: iq
    integer(IK) :: iw
    integer(IK) :: rest
    integer(IK) :: xtry
    integer(IK) :: cp
    integer(IK) :: wp
    integer(IK) :: cq
    integer(IK) :: wq
  end type mol_block_
!
  interface
    include 'dgemm.h'
  end interface
!
  integer(IK), parameter :: DEF_maxiter   = 1000_IK
  real(RK), parameter    :: DEF_threshold = 1.0E-12_RK
!
contains
!
  pure elemental function b2b_(d, x, y, p, q, w, b) result(res)
    integer(IK), intent(in)     :: d
!!! spatial dimension
    integer(IK), intent(in)     :: x, y, p, q, w
!!! pointer to x, y, p, q, w
    type(mol_block), intent(in) :: b
    type(mol_block_)            :: res
    res%m = b%m
    res%n = b%n
    res%g = b%g
    res%f = b%f
    res%mg = res%m * res%g
    res%mn = res%m * res%n
    res%df = d * res%f
    res%dm = d * res%m
    res%ff = res%f * res%f
    res%gg = res%g * res%g
    res%ffg = res%ff * res%g
    res%dmn = d * res%mn
    res%dmg = d * res%mg
    res%m_f = res%m - res%f
    res%ffgg = res%ff * res%gg
    res%proc_f = Procrustes_worksize(res%f)
    res%blksiz = res%ff + res%proc_f
    res%dxm_dxf = d * res%m_f
    res%ix = x
    res%iy = y
!!! rotation matrix P(g,g)
    res%ip = p
!!! rotation matrix Q(f,f,g)
    res%iq = q
!!! covariance matrices and residue matrix, REST(g, g), XTRY(f, f, g, g)
    res%iw = w
    res%rest = 1
    res%xtry = res%rest + res%gg
!!! covariance matrix for P estimation, CP(g, g), w(*)
    res%cp = res%xtry + res%ffgg
    res%wp = res%cp + res%gg
!!! covariance matrix for Q estimation, CQ(f, f, g), w(*)
    res%cq = res%xtry + res%ffgg
    res%wq = res%cq + res%ff
  end function b2b_
!
!| Calculate work array size for d*d matrix.
  pure function block_lower_bound_worksize(d, m, n, free_m_indices, free_n_indices) result(res)
    integer(IK), intent(in)           :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)           :: m
    !! matrix dimension 2.
    integer(IK), intent(in)           :: n
    !! matrix dimension 3.
    integer(IK), intent(in)           :: free_m_indices(:)
    !! rotation index, must be SIZE(free_indices) <= m
    integer(IK), intent(in)           :: free_n_indices(:)
    !! rotation index, must be SIZE(free_indices) <= n
    integer(IK)                       :: res
    integer(IK)                       :: f, g
!
    if (d < 1 .or. m < 1 .or. n < 1) then
      res = 0
    else
      f = SIZE(free_m_indices)
      g = SIZE(free_n_indices)
      res = 2 + 3 * d * m * n &
     &    + MAX(d * d * 2 + Kabsch_worksize(d), &
     &          (m * n)**2 + g * (2 * f * f + procrustes_worksize(f)))
    end if
!
  end function block_lower_bound_worksize
!
!| Calculate min_{P,Q} |X-PYQ|^2, Q is block rotation matrix
  pure subroutine block_lower_bound(d, s, b, X, Y, w, maxiter, threshold, R0)
    integer(IK), intent(in)           :: d
    !! matrix dimension 1.
    integer(IK), intent(in)           :: s
    !! matrix dimension 1.
    type(mol_block), intent(in)       :: b(s)
!   integer(IK), intent(in)           :: m
!   !! matrix dimension 2.
!   integer(IK), intent(in)           :: n
!   !! matrix dimension 3.
!   integer(IK), intent(in)           :: free_m_indices(:)
!   !! rotation index, must be SIZE(free_indices) <= m
!   integer(IK), intent(in)           :: free_n_indices(:)
    !! rotation index, must be SIZE(free_indices) <= n
    real(RK), intent(in)              :: X(*)
    !! reference d*sum_i (mi*ni) array
    real(RK), intent(in)              :: Y(*)
    !! target d*sum_i (mi*ni) array
    real(RK), intent(inout)           :: w(*)
    !! work array, must be larger than lower_bound_worksize(d, n, free_indices)
    integer(IK), intent(in), optional :: maxiter
    !! iteration limit, default = 1000
    real(RK), intent(in), optional    :: threshold
    !! iteration limit, default = 1E-12
    real(RK), intent(in), optional    :: R0(*)
    !! initial Rotation matrix
    integer(IK)                       :: maxiter_
    real(RK)                          :: threshold_
    real(RK)                          :: conv_(s), conv, prev
    type(mol_block_)                  :: b_(s)
    integer(IK)                       :: dd
    integer(IK)                       :: cf, r, p, q
    integer(IK)                       :: ix, iy
    integer(IK)                       :: i, j, iw
!
    if (d < 1 .or. s < 1) then
      w(1) = rmsd(d, SUM(b%m * b%n), X, Y)
      return
    end if
!
    dd = d * d
!   mg = m * g
!   mn = m * n
!   dm = d * m
!   df = d * f
!   ff = f * f
!   gg = g * g
!   dmn = d * m * n
!
!   ix   = 1          ! copy of X
!   iy   = ix   + dmn ! copy of Y
!   cf   = iy   + dmn ! covariance matrix X^TYPQ
!   r    = cf   + dd  ! rotation matrix R
!
!   p    = r    + dd     ! rotation matrix R
!   q    = p    + gg     ! rotation matrix R
!   iw   = q    + ff * g
    block
      integer(IK) :: iix, iiy
      iix = ix
      iiy = iy
      do i = 1, s
        !b_(i) = b2b_(d, b(i))
      end do
    end block
!
!   call get_C_fix(d, m, n, g, w(ix), w(iy), w(cf))
!
    maxiter_   = optarg(maxiter,   DEF_maxiter)
    threshold_ = optarg(threshold, DEF_threshold) * SUM(b_%dmn)
!
!   call calc_swap_indices(m, f, free_m_indices, m_indices)
!   call calc_swap_indices(n, g, free_n_indices, n_indices)
!   call pack_matrix(d, m, n, m_indices, n_indices, X, w(ix))
!   call pack_matrix(d, m, n, m_indices, n_indices, Y, w(iy))
!   call copy(dmn, w(iy), w(iz))
!
    conv_ = -RHUGE
!
    if (PRESENT(R0)) then
      call copy(dd, R0, w(r))
    else
      call eye(d, w(r))
    end if
!
    do i = 1, maxiter_
      prev = conv
      do concurrent(j=1:s)
     !  call update_PQ(d, b_(j), w(b_(j)%ix), w(b_(j)%iy), w(r), maxiter_, threshold_, &
     ! &               w(b_(j)%ip), w(b_(j)%iq), conv_(j), W(b_(j)%iw))
      end do
      !call update_R(d, s, b_, w(ix), w(iy), w(cf), w(r), W(iw))
      conv = SUM(conv_)
      if (ABS(conv - prev) < threshold_) exit
    end do
!
!   call R_rotation(d, mn, w(r), w(iy), w(cf))
!   call PQ_rotation(d, b_, w(p), w(q), w(iy), w(cf))
!   w(1) = rmsd(d, mn, w(ix), w(iy))
!
  contains
!
    pure subroutine pack_matrix(d, m, n, m_indices, n_indices, source, dest)
      integer(IK), intent(in) :: d, m, n, m_indices(m), n_indices(n)
      real(RK), intent(in)    :: source(d, m, n)
      real(RK), intent(inout) :: dest(d, m, n)
      integer(IK)             :: i, j, k
!
      do concurrent(k=1:n, j=1:m, i=1:d)
        dest(i, j, k) = source(i, m_indices(j), n_indices(k))
      end do
!
    end subroutine pack_matrix
!
    pure subroutine calc_swap_indices(n, f, free_indices, res)
      integer(IK), intent(in)    :: n, f, free_indices(f)
      integer(IK), intent(inout) :: res(n)
      integer(IK)                :: i, fi, n_f
      do concurrent(i=1:f)
        res(i) = free_indices(i)
      end do
      n_f = n - f
      fi  = f + 1
      do i = 1, n
        if (ANY(i == res(:f))) cycle
        res(fi) = i
        if (fi == n) exit
        fi = fi + 1
      end do
    end subroutine calc_swap_indices
!
  end subroutine block_lower_bound
!
  pure subroutine get_C_fix(d, m, n, g, X, Y, CF)
    integer(IK), intent(in) :: d, m, n, g
    real(RK), intent(in)    :: X(d, m, n), Y(d, m, n)
    real(RK), intent(inout) :: CF(d, d)
    integer(IK)             :: mn_g
      mn_g = m * (n - g)
      if (mn_g < 1) then
        call zfill(d * d, CF) ; return
      end if
      call DGEMM('N', 'T', d, d, mn_g, ONE, X(:, :, g + 1), d, Y(:, :, g + 1), d, ZERO, CF, d)
  end subroutine get_C_fix
!
  pure subroutine update_R(d, s, b, X, Y, CF, R, Z, W)
    integer(IK), intent(in)      :: d, s
    type(mol_block_), intent(in) :: b(s)
    real(RK), intent(in)         :: X(*), Y(*)
    !real(RK), intent(in)         :: P(g, g), Q(f, f, g), CF(d, d)
    real(RK), intent(in)         :: CF(d, d)
    real(RK), intent(inout)      :: R(d, d), Z(*), W(*)
    integer(IK)                  :: i
!!    R rotation
!!    M = X(d,mn)@Z^T(mn,d)@P^T@Q^T
      do concurrent(i=1:s)
        call copy(b(i)%dmg, Y(b(i)%ix), Z(b(i)%ix))
        call PQ_rotation(d, b(i), w(b(i)%ip), w(b(i)%iq), Z(b(i)%ix), W(b(i)%iw))
      end do
!!    W = CF + YPQXT
      call copy(d * d, CF, W)
      call DGEMM('N', 'T', d, d, SUM(b%mg), ONE, X, d, Z, d, ONE, W, d)
!!    get optimal R
      call Kabsch(d, W, R, W(d * d + 1))
  end subroutine update_R
!
  pure subroutine update_XRYT(d, b, X, Y, R, XTRY, RESTR)
    !integer(IK), intent(in) :: d, m, n, f, g
    integer(IK), intent(in)      :: d
    type(mol_block_), intent(in) :: b
    real(RK), intent(in)         :: X(d, b%m, b%n), Y(d, b%m, b%n), R(d, d)
    real(RK), intent(inout)      :: XTRY(b%f, b%f, b%g, b%g), RESTR(b%g, b%g)
    integer(IK)                  :: i, j
!
    if (b%f > 0) then
      do concurrent(i=1:b%g, j=1:b%g)
        block
          real(RK) :: temp(b%df)
          call DGEMM('T', 'N', b%f, d, d, ONE, X(1, 1, i), d, R, d, ZERO, temp, b%f)
          call DGEMM('N', 'N', b%f, b%f, d, ONE, temp, b%f, Y(1, 1, j), d, ZERO, XTRY(1, 1, i, j), b%f)
        end block
      end do
    end if
!
    if (b%m_f < 1) return
!
    do concurrent(i=1:b%g, j=1:b%g)
      block
        real(RK) :: RY(b%dxm_dxf)
        call DGEMM('N', 'N', d, b%m_f, d, ONE, R, d, Y(1, b%f + 1, j), d, ZERO, RY, d)
        RESTR(i, j) = ddot(b%dxm_dxf, X(1, b%f + 1, i), RY)
      end block
    end do
!
  end subroutine update_XRYT
!
  pure subroutine update_PQ(d, b, X, Y, R, Q, P, maxiter, threshold, trace, W)
    !integer(IK), intent(in) :: d, m, n, f, g, maxiter
    integer(IK), intent(in)      :: d, maxiter
    type(mol_block_), intent(in) :: b
    real(RK), intent(in)         :: X(*), Y(*), R(d, d), threshold
    real(RK), intent(inout)      :: P(b%g, b%g), Q(b%f, b%f, b%g)
    real(RK), intent(inout)      :: trace, W(*)
    real(RK)                     :: prev
    !integer(IK)                  :: rest, xtry, cp, cq, wp, wq
    integer(IK)                  :: i, j, l
!
!!! covariance matrices and residue matrix, REST(g, g), XTRY(f, f, g, g)
!   rest = 1
!   xtry = rest + b%gg
!
!!! covariance matrix for P estimation, CP(g, g), w(*)
!   cp = xtry + b%ffgg
!   wp = cp + b%gg
!
!!! covariance matrix for Q estimation, CQ(f, f, g), w(*)
!   cq = xtry + b%ffgg
!   wq = cq + b%ff
!
    trace = -RHUGE
!
    call update_XRYT(d, b, X, Y, R, w(b%xtry), w(b%rest))
!
    call eye(b%g, P)
    do concurrent(i=1:b%g)
      call eye(b%f, Q(1, 1, i))
    enddo
!
    do l = 1, maxiter
!
      prev = trace
!
!!!   CP = RESTR
      call copy(b%gg, w(b%rest), w(b%cp))
      do concurrent(i=1:b%g, j=1:b%g)
!!!     CP(i,j) = trace(Qi@Cij) = ddot(Qi^T, Cij)
        block
          integer(IK) :: icp, ixtry
          icp = b%cp + i + (j - 1) * b%g - 1
          ixtry = b%xtry + b%ff * ((i - 1) + b%g * (j - 1))
          w(icp) = w(icp) + ddot(b%ff, Q(1, 1, i), w(ixtry))
        end block
      end do
!!!   Get P = UV^T
      call Procrustes(b%g, w(b%cp), P, w(b%wp))
!!!   tr(P^T, CP) = ddot(P, CP)
      trace = ddot(b%gg, P, w(b%cp))
      if (ABS(trace - prev) < threshold) exit
!
      do concurrent(j=1:b%g)
        block
          integer(IK) :: icq, iwq, ixtry
          icq = b%cq + b%blksiz * (j - 1)
          iwq = b%wq + b%blksiz * (j - 1)
          ixtry = b%xtry + b%ff * (j - 1)
!!        CQ^T = sum_i Pi * Cji
          call zfill(b%ff, w(icq))
          do i = 1, b%g
            call add(b%ff, P(j, i), w(ixtry), w(icq))
            ixtry = ixtry + b%ffg
          end do
!!!       Get Qi^T = UV^T
          call Procrustes(b%f, w(icq), Q(1, 1, j), w(iwq))
        end block
      end do
    end do
!
  end subroutine update_PQ
!
  pure subroutine R_rotation(d, n, R, X, W)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: R(*)
    real(RK), intent(inout) :: X(*), W(*)
    call DGEMM('N', 'N', d, n, d, ONE, R, d, X, d, ZERO, W, d)
    call copy(d * n, W, X)
  end subroutine R_rotation
!
  pure subroutine PQ_rotation(d, b, P, Q, X, W)
    integer(IK), intent(in)      :: d
    type(mol_block_), intent(in) :: b
    real(RK), intent(in)         :: P(*), Q(*)
    real(RK), intent(inout)      :: X(*), W(*)
    integer(IK)                  :: i, ff, dm, df, dmg, dm_f
!
    call DGEMM('N', 'N', b%dm, b%g, b%g, ONE, X, b%dm, P, b%g, ZERO, W, b%dm)
!
    do concurrent(i=1:b%g)
      block
        integer(IK) :: iq, ix
        ix = (i - 1) * b%dm + 1
        iq = (i - 1) * b%ff + 1
        call DGEMM('N', 'T', d, b%f, b%f, ONE, W(ix), d, q(iq), b%f, ZERO, X(ix), d)
        ix = ix + df
        call copy(b%dxm_dxf, W(ix), X(ix))
      end block
    end do
!
  end subroutine PQ_rotation
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
  pure subroutine add(d, c, source, dest)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: c, source(*)
    real(RK), intent(inout) :: dest(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dest(i) = dest(i) + c * source(i)
    end do
  end subroutine add
!
  pure function ddot(d, X, Y) result(res)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: X(*), Y(*)
    real(RK)                :: res
    integer(IK)             :: i
    res = ZERO
    do i = 1, d
      res = res + X(i) * Y(i)
    end do
  end function ddot
!
  pure subroutine eye(d, X)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: X(d, *)
    integer(IK)             :: i, j
    do concurrent(i=1:d, j=1:d)
      X(i, j) = MERGE(ONE, ZERO, i == j)
    end do
  end subroutine eye
!
  pure subroutine zfill(d, dest)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: dest(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dest(i) = ZERO
    end do
  end subroutine zfill
!
end module mod_block_lower_bound
