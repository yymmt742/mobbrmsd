!| Calculate min_{R,P,Q1,...,Q_n} |X-R@Y@f(P)@g(Q1,Q2,...,Qn)|^2 for R\in\R^{d,d}, P\in \R^{n,n} and Qi\in \R^{m,m}.
!  f(P) = I x P ( I is the identity and \in R^{m,m} ).
!  g(Q1,Q2,...,Qn) = Diag(Q1 Q2 ... Qn).
!  Here, R@R^T=I, det(R)=1, P@P^T=I, and Q@Q^T=I are satisfied.
module mod_block_lower_bound
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_optarg
  use mod_rmsd
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
  subroutine block_lower_bound(d, m, n, free_m_indices, free_n_indices, X, Y, w, maxiter, threshold)
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
    integer(IK)                       :: maxiter_
    real(RK)                          :: threshold_
!   real(RK)                          :: cov(SIZE(free_n_indices) * m, SIZE(free_n_indices) * m)
!   real(RK)                          :: mrot(SIZE(free_m_indices), SIZE(free_m_indices), SIZE(free_n_indices))
!   real(RK)                          :: nrot(SIZE(free_n_indices), SIZE(free_n_indices))
    integer(IK)                       :: n_indices(n), m_indices(m)
    integer(IK)                       :: f, g, mg, mn, dd, dm, df, ff, dfg, dmg, dmn, mmn, mmgg
    integer(IK)                       :: conv, prev
    integer(IK)                       :: ix, iy, iz, rcov, rrot, cov, mrot, nrot
    integer(IK)                       :: i, j, iw1, iw2
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
    dmn = d * m * n
    mmn = m * m * n
    mmgg = mg * mg
    dfg = d * f * g
    dmg = d * m * g
!
    conv = 1
    prev = conv + 1
    ix   = prev + 1      ! copy of X
    iy   = ix   + dmn    ! copy of Y
    iz   = iy   + dmn    ! Z = P@Y
    rcov = iz   + dmn    ! covariance matrix X^TYPQ
    rrot = rcov + dd     ! rotation matrix R
    iw1  = rrot + dd
    cov  = iz   + dmn    ! covariance matrix XRY^T
    mrot = cov  + mmgg   ! rotation matrix P
    nrot = mrot + ff * g ! rotation matrix Q
    iw2 = nrot + g * g
!
    maxiter_   = optarg(maxiter,   DEF_maxiter)
    threshold_ = optarg(threshold, DEF_threshold) * dmn
!
    call calc_swap_indices(m, f, free_m_indices, m_indices)
    call calc_swap_indices(n, g, free_n_indices, n_indices)
!
    call pack_matrix(d, m, n, m_indices, n_indices, X, w(ix))
    call pack_matrix(d, m, n, m_indices, n_indices, Y, w(iy))
!
    do concurrent(j=0:dmn-1)
      w(iz + j) = w(iy + j)
    end do
!
    w(conv) = -RHUGE
!
    do i = 1, maxiter_
!
      w(prev) = w(conv)
!
!!    M = X^T@Z
      call DGEMM('T', 'N', mg, mg, d, ONE, w(ix), d, w(iz), d, ZERO, w(cov), mg)
      call block_rotation(f, g, m, mg, w(cov), w(mrot), w(nrot), maxiter_, threshold_)
!
!!    R rotation
!!    M = X(d,mn)@Z^T(mn,d)@P^T@Q^T
      w(iz:iz + dmg - 1) = w(iy:iy + dmg - 1)
      call P_rotation(d, m, n, f, g, w(mrot), w(iz))
      call Q_rotation(dm, n, g, w(nrot), w(iz))
      call DGEMM('N', 'T', d, d, mn, ONE, w(ix), d, w(iz), d, ZERO, w(rcov), d)
!!    get optimal R
      call Kabsch(d, w(rcov), w(rrot), w(iw1))
!
!!    trace(MR) = trace(MR^T)
!
      w(conv) = ZERO
      do j = 0, dd - 1
        w(conv) = w(conv) + w(rcov + j) * w(rrot + j)
      end do
!
      if (ABS(w(prev) - w(conv)) < threshold_) exit
!
!!    Z = R@Y
      call DGEMM('N', 'N', d, mn, d, ONE, w(rrot), d, w(iy), d, ZERO, w(iz), d)
!
    enddo
!
    w(iz:iz + dmg - 1) = w(iy:iy + dmg - 1)
    call P_rotation(d, m, n, f, g, w(mrot), w(iz))
    call Q_rotation(dm, n, g, w(nrot), w(iz))
    w(conv) = rmsd(d, mn, w(ix), w(iz))
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
  end subroutine block_lower_bound
!
  pure subroutine block_rotation(f, g, m, mg, cov, mrot, nrot, maxiter, threshold)
    integer(IK), intent(in) :: f, g, m, mg, maxiter
    real(RK), intent(in)    :: cov(mg, mg)
    real(RK), intent(in)    :: threshold
    real(RK), intent(inout) :: mrot(f, f, g), nrot(g, g)
    real(RK)                :: mcov(f, f), ncov(g, g), ncov0(g, g)
    real(RK)                :: w(procrustes_worksize(f))
    real(RK)                :: prev, conv
    integer(IK)             :: i, j, k
!
    do concurrent(j=1:g, k=1:g)
      block
        integer(IK) :: jm, km
        jm = (j - 1) * m
        km = (k - 1) * m
        ncov0(j, k) = SUM(cov(jm + f + 1:jm + m, km + 1:km + f))&
       &            + SUM(cov(jm + 1:jm + m, km + f + 1:km + m))
      end block
    end do
!
    prev = -RHUGE
    conv = -RHUGE
!
    do concurrent(i=1:g, j=1:g)
      nrot(i, j) = MERGE(1, 0, i == j)
    end do
!
    do concurrent(i=1:f, j=1:f, k=1:g)
      mrot(i, j, k) = MERGE(1, 0, i == j)
    end do
!
    do i = 1, maxiter
!
      if (prev < conv) prev = conv
!
!!    P rotation
!
      do concurrent(j=1:g)
        mcov = ZERO
        do k = 1, g
          mcov(:, :) = mcov(:, :) &
         &           + nrot(k, j) * cov((k - 1) * m + 1:(k - 1) * m + f, (k - 1) * m + 1:(k - 1) * m + f)
        end do
        call Procrustes(f, mcov, mrot(:, :, j), w)
      end do
!
!!    Q rotation
      do concurrent(j=1:g, k=1:g)
        block
          integer(IK) :: jm, km
          jm = (j - 1) * m
          km = (k - 1) * m
          ncov(j, k) = SUM(MATMUL(cov(jm + 1:jm + f, km + 1:km + f), TRANSPOSE(mrot(:, :, k))))
        end block
      end do
!
      call Procrustes(g, ncov, nrot, w)
!
      conv = SUM(ncov * nrot)
      if (ABS(prev - conv) < threshold) exit
!
    end do
!
  end subroutine block_rotation
!
  pure subroutine P_rotation(d, m, n, f, g, r, z)
    integer(IK), intent(in) :: d, m, n, f, g
    real(RK), intent(in)    :: r(f, f, g)
    real(RK), intent(inout) :: z(d, m, n)
    real(RK)                :: t(d, f)
    integer(IK)             :: i
!
    do i = 1, g
      t = MATMUL(z(:, :f, i), TRANSPOSE(r(:, :, i)))
      z(:, :f, i) = t
      !z(:, :f, i) = MATMUL(z(:, :f, i), TRANSPOSE(r(:, :, i)))
    end do
!
  end subroutine P_rotation
!
  pure subroutine Q_rotation(dm, n, g, r, z)
    integer(IK), intent(in) :: dm, n, g
    real(RK), intent(in)    :: r(g, g)
    real(RK), intent(inout) :: z(dm, n)
!
    z(:, :g) = MATMUL(z(:, :g), TRANSPOSE(r))
!
  end subroutine Q_rotation
!
  pure subroutine calc_swap_indices(n, f, free_indices, res)
    integer(IK), intent(in)    :: n, f, free_indices(f)
    integer(IK), intent(inout) :: res(n)
    integer(IK)                :: i, j, g
!
    do concurrent(i=1:f)
      res(i) = free_indices(i)
    end do
!
    g = n - f
    j = f + 1
!
    do i = 1, n
      if (ANY(i == res(:f))) cycle
      res(j) = i
      if (j == n) exit
      j = j + 1
    end do
!
  end subroutine calc_swap_indices
!
  pure function objective(d, n, m, f, g, P, Q, R) result(res)
    integer(IK), intent(in) :: d, n, m, f, g
    real(RK), intent(in)    :: P(d, d), Q(g, g), R(f, f, g)
    real(RK), intent(in)    :: x(*), y(*)
    real(RK)                :: YTQRX()
    real(RK)                :: res
!
  end function objective
!
end module mod_block_lower_bound
