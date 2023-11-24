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
  interface
    include 'dgemm.h'
  end interface
!
  integer(IK), parameter :: DEF_maxiter   = 1000_IK
  real(RK), parameter    :: DEF_threshold = 1.0E-12_RK
!
contains
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
      res = 2 + 3 * d * m * n + MAX(d * d * 2 + Kabsch_worksize(d), (m * n)**2 + g * (2 * f * f + procrustes_worksize(f)))
    end if
!
  end function block_lower_bound_worksize
!
!| Calculate min_{P,Q} |X-PYQ|^2, Q is block rotation matrix
  subroutine block_lower_bound(d, m, n, free_m_indices, free_n_indices, X, Y, w, maxiter, threshold, R0)
    integer(IK), intent(in)           :: d
    !! matrix dimension 1.
    integer(IK), intent(in)           :: m
    !! matrix dimension 2.
    integer(IK), intent(in)           :: n
    !! matrix dimension 3.
    integer(IK), intent(in)           :: free_m_indices(:)
    !! rotation index, must be SIZE(free_indices) <= m
    integer(IK), intent(in)           :: free_n_indices(:)
    !! rotation index, must be SIZE(free_indices) <= n
    real(RK), intent(in)              :: X(*)
    !! reference d*m*n array
    real(RK), intent(in)              :: Y(*)
    !! target d*m*n array
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
    integer(IK)                       :: n_indices(n), m_indices(m)
    integer(IK)                       :: f, g, mg, mn, dd, dm, df, ff, gg, dmn
    integer(IK)                       :: cf, r, p, q
    integer(IK)                       :: ix, iy, iz
    integer(IK)                       :: conv, prev
    integer(IK)                       :: i, iw
!
    f = SIZE(free_m_indices)
    g = SIZE(free_n_indices)
!
    if (f < 1 .or. g < 1) then
      w(1) = rmsd(d, m * n, X, Y)
      return
    end if
!
    mg = m * g
    mn = m * n
    dm = d * m
    dd = d * d
    df = d * f
    ff = f * f
    gg = g * g
    dmn = d * m * n
!
    conv = 1
    prev = conv + 1
    ix   = prev + 1      ! copy of X
    iy   = ix   + dmn    ! copy of Y
    iz   = iy   + dmn    ! Z = P@Y
    cf   = iz   + dmn    ! covariance matrix X^TYPQ
    r    = cf   + dd     ! rotation matrix R
    p    = r    + dd     ! rotation matrix R
    q    = p    + gg     ! rotation matrix R
    iw   = q    + ff * g
!
    maxiter_   = optarg(maxiter,   DEF_maxiter)
    threshold_ = optarg(threshold, DEF_threshold) * dmn
!
    call calc_swap_indices(m, f, free_m_indices, m_indices)
    call calc_swap_indices(n, g, free_n_indices, n_indices)
    call pack_matrix(d, m, n, m_indices, n_indices, X, w(ix))
    call pack_matrix(d, m, n, m_indices, n_indices, Y, w(iy))
    call copy(dmn, w(iy), w(iz))
!
    w(conv) = -RHUGE
!
    call get_C_fix(d, m, n, g, w(ix), w(iy), w(cf))
!
    if (PRESENT(R0)) then
      call copy(dd, R0, w(r))
    else
      call eye(d, w(r))
    end if
!
    do i = 1, maxiter_
      w(prev) = w(conv)
      call update_PQ(d, m, n, f, g, w(ix), w(iy), w(r), maxiter_, threshold_, w(p), w(q), w(conv), W(iw))
      call update_R(d, m, n, f, g, w(ix), w(iy), w(p), w(q), w(cf), w(r), W(iw))
      if (ABS(w(conv) - w(prev)) < threshold_) exit
    end do
!
    call R_rotation(d, mn, w(r), w(iy), w(iz))
    call PQ_rotation(d, m, n, f, g, w(p), w(q), w(iy), w(iz))
    w(conv) = rmsd(d, mn, w(ix), w(iy))
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
  pure subroutine update_R(d, m, n, f, g, X, Y, P, Q, CF, R, W)
    integer(IK), intent(in) :: d, m, n, f, g
    real(RK), intent(in)    :: X(d, m, n), Y(d, m, n)
    real(RK), intent(in)    :: P(g, g), Q(f, f, g), CF(d, d)
    real(RK), intent(inout) :: R(d, d), W(*)
    integer(IK)             :: dd, mg, dmg, iz, ic, iw
!!    R rotation
!!    M = X(d,mn)@Z^T(mn,d)@P^T@Q^T
      dd = d * d
      mg = m * g
      dmg = d * mg
      ic = 1
      iz = ic + dd
      iw = iz + dmg
      call copy(dmg, Y, W(iz))
      call PQ_rotation(d, m, n, f, g, P, Q, W(iz), W(iw))
      call copy(dd, CF, W(ic))
      call DGEMM('N', 'T', d, d, mg, ONE, X, d, W(iz), d, ONE, W(ic), d)
!!    get optimal R
      call Kabsch(d, W(ic), R, W(iz))
  end subroutine update_R
!
  pure subroutine update_XRYT(d, m, n, f, g, X, Y, R, XTRY, RESTR)
    integer(IK), intent(in) :: d, m, n, f, g
    real(RK), intent(in)    :: X(d, m, n), Y(d, m, n), R(d, d)
    real(RK), intent(inout) :: XTRY(f, f, g, g), RESTR(g, g)
    integer(IK)             :: i, j, m_f
!
    if (f > 0) then
      do concurrent(i=1:g, j=1:g)
        block
          real(RK) :: temp(d * f)
          call DGEMM('T', 'N', f, d, d, ONE, X(1, 1, i), d, R, d, ZERO, temp, f)
          call DGEMM('N', 'N', f, f, d, ONE, temp, f, Y(1, 1, j), d, ZERO, XTRY(1, 1, i, j), f)
        end block
      end do
    end if
!
    m_f = m - f; if (m_f < 1) return
!
    do concurrent(i=1:g, j=1:g)
      block
        real(RK) :: RY(d * m_f)
        call DGEMM('N', 'N', d, m_f, d, ONE, R, d, Y(1, f + 1, j), d, ZERO, RY, d)
        RESTR(i, j) = ddot(d * m_f, X(1, f + 1, i), RY)
      end block
    end do
!
  end subroutine update_XRYT
!
  pure subroutine update_PQ(d, m, n, f, g, X, Y, R, maxiter, threshold, P, Q, trace, W)
    integer(IK), intent(in) :: d, m, n, f, g, maxiter
    real(RK), intent(in)    :: X(d, m, n), Y(d, m, n), R(d, d), threshold
    real(RK), intent(inout) :: P(g, g), Q(f, f, g), trace, W(*)
    real(RK)                :: prev
    integer(IK)             :: rest, xtry, cp, cq, wp, wq, bs, gg, ff
    integer(IK)             :: i, j, l
!
    gg = g * g
    ff = f * f
!
!!! covariance matrices and residue matrix, REST(g, g), XTRY(f, f, g, g)
    rest = 1
    xtry = rest + gg
!
!!! covariance matrix for P estimation, CP(g, g), w(*)
    cp = xtry + ff * gg
    wp = cp + gg
!
!!! covariance matrix for Q estimation, CQ(f, f, g), w(*)
    cq = xtry + ff * gg
    wq = cq + ff
    bs = ff + Procrustes_worksize(f)
!
    trace = -RHUGE
!
    call update_XRYT(d, m, n, f, g, X, Y, R, w(xtry), w(rest))
!
    call eye(g, p)
    do concurrent(i=1:g)
      call eye(f, Q(1, 1, i))
    enddo
!
    do l = 1, maxiter
!
      prev = trace
!
!!!   CP = RESTR
      call copy(gg, w(rest), w(cp))
      do concurrent(i=1:g, j=1:g)
!!!     CP(i,j) = trace(Qi@Cij) = ddot(Qi^T, Cij)
        block
          integer(IK) :: icp, ixtry
          icp = cp + i + (j - 1) * g - 1
          ixtry = xtry + ff * ((i - 1) + g * (j - 1))
          w(icp) = w(icp) + ddot(ff, Q(1, 1, i), w(ixtry))
        end block
      end do
!!!   Get P = UV^T
      call Procrustes(g, w(cp), P, w(wp))
!!!   tr(P^T, CP) = ddot(P, CP)
      trace = ddot(gg, P, w(cp))
      if (ABS(trace - prev) < threshold) exit
!
      do concurrent(j=1:g)
        block
          integer(IK) :: icq, iwq, ixtry
          icq = cq + bs * (j - 1)
          iwq = wq + bs * (j - 1)
          ixtry = xtry + ff * (j - 1)
!!        CQ^T = sum_i Pi * Cji
          call zfill(ff, w(icq))
          do i = 1, g
            call add(ff, P(j, i), w(ixtry), w(icq))
            ixtry = ixtry + g * ff
          end do
!!!       Get Qi^T = UV^T
          call Procrustes(f, w(icq), Q(1, 1, j), w(iwq))
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
  pure subroutine PQ_rotation(d, m, n, f, g, P, Q, X, W)
    integer(IK), intent(in) :: d, m, n, f, g
    real(RK), intent(in)    :: P(*), Q(*)
    real(RK), intent(inout) :: X(*), W(*)
    integer(IK)             :: i, ff, dm, df, dmg, dm_f
!
    ff = f * f
    dm = d * m
    df = d * f
    dmg = dm * g
    dm_f = d * (m - f)
!
    call DGEMM('N', 'N', dm, g, g, ONE, X, dm, P, g, ZERO, W, dm)
!
    do concurrent(i=1:g)
      block
        integer(IK) :: iq, ix
        ix = (i - 1) * dm + 1
        iq = (i - 1) * ff + 1
        call DGEMM('N', 'T', d, f, f, ONE, W(ix), d, q(iq), f, ZERO, X(ix), d)
        ix = ix + df
        call copy(dm_f, W(ix), X(ix))
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
