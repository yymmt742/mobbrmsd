!| Calculate min_{R,P} |X-RYP|^2
!  Here, RR^T=PP^T=I and det(R)=1 are satisfied.
module mod_lower_bound
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_optarg
  use mod_rmsd
  use mod_cov
  use mod_Kabsch
  use mod_Procrustes
  implicit none
  private
  public :: lower_bound_worksize, lower_bound, &
          & block_lower_bound_worksize, block_lower_bound
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
  pure function lower_bound_worksize(d, n, free_indices) result(res)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)       :: n
    !! matrix row dimension.
    integer(IK), intent(in)       :: free_indices(:)
    !! rotation index, must be SIZE(free_indices) <= n
    integer(IK)                   :: res
    integer(IK)                   :: f
!
    if (d < 1 .or. n < 1) then
      res = 0
    else
      f = SIZE(free_indices)
      res = 2 + 3 * d * n + MAX(d * d * 2 + Kabsch_worksize(d), f * f * 2 + procrustes_worksize(f))
    end if
!
  end function lower_bound_worksize
!
!| Calculate min_{P,Q} |X-PYQ|^2
  pure subroutine lower_bound(d, n, free_indices, X, Y, w, maxiter, threshold)
  !pure subroutine lower_bound(d, n, free_indices, X, Y, w, maxiter, threshold)
    integer(IK), intent(in)           :: d
    !! matrix collumn dimension.
    integer(IK), intent(in)           :: n
    !! matrix row dimension.
    integer(IK), intent(in)           :: free_indices(:)
    !! rotation index, must be SIZE(free_indices) <= n
    real(RK), intent(in)              :: X(*)
    !! reference d*n array
    real(RK), intent(in)              :: Y(*)
    !! target d*n array
    real(RK), intent(inout)           :: w(*)
    !! work array, must be larger than lower_bound_worksize(d, n, free_indices)
    integer(IK), intent(in), optional :: maxiter
    !! iteration limit, default = 1000
    real(RK), intent(in), optional    :: threshold
    !! iteration limit, default = 1E-12
    integer(IK)                       :: i, j, maxiter_
    real(RK)                          :: threshold_
    integer(IK)                       :: f, dn, dd, df, ff
    integer(IK)                       :: conv, prev
    integer(IK)                       :: ix, iy, iz, drot, ncov, nrot
    integer(IK)                       :: iw1, iw2
!
    f = SIZE(free_indices)
    if(f<1) return
!
    dn = d * n
    dd = d * d
    df = d * f
    ff = f * f
!
    conv = 1
    prev = conv + 1
    ix   = prev + 1  ! copy of X
    iy   = ix   + dn ! copy of Y
    iz   = iy   + dn ! Z = P@Y
    drot = iz   + dn ! rotation matrix P
    iw1  = drot + dd
    nrot = iz   + dn ! rotation matrix Q
    ncov = nrot + ff ! covariance matrix, X^TY
    iw2  = nrot + ff
!
    maxiter_   = optarg(maxiter,   DEF_maxiter)
    threshold_ = optarg(threshold, DEF_threshold) ** 2 * n
!
    call pack_matrix(d, n, f, free_indices, X, w(ix))
    call pack_matrix(d, n, f, free_indices, Y, w(iy))
!
    do concurrent(j=0:dn - 1)
      w(iz + j) = w(iy + j)
    end do
!
    w(conv) = sd(d, n, w(ix), w(iy))
!
    do i = 1, maxiter_
!
      w(prev) = w(conv)
!
      call R_rotation(d, n, w(ix), w(iy), w(iz), w(drot), w(iw1))
!
!!    M = X^T@Z
      call DGEMM('T', 'N', f, f, d, ONE, w(ix), d, w(iz), d, ZERO, w(ncov), f)
!!    calc P from M
      call Procrustes(f, w(ncov), w(nrot), w(iw2))
!!    Z' = Z@P
      call DGEMM('N', 'T', d, f, f, ONE, w(iz), d, w(nrot), f, ZERO, w(iw2), d)
!!    Z = Z'
      do concurrent(j=0:df - 1)
        w(iz + j) = w(iw2 + j)
      end do
!
      w(conv) = sd(d, n, w(ix), w(iz))
      if (ABS(w(prev) - w(conv)) < threshold_) exit
!
!!    Z = Y
      do concurrent(j=0:dn-1)
        w(iz + j) = w(iy + j)
      enddo
!!    Z' = Z@P
      call DGEMM('N', 'T', d, f, f, ONE, w(iz), d, w(nrot), f, ZERO, w(iw2), d)
!!    Z = Z'
      do concurrent(j=0:df - 1)
        w(iz + j) = w(iw2 + j)
      end do
!
    enddo
!
    w(conv) = SQRT(w(conv) / n)
!
  contains
!
    pure subroutine pack_matrix(d, n, f, free_indices, source, dest)
      integer(IK), intent(in) :: d, n, f, free_indices(f)
      real(RK), intent(in)    :: source(d, n)
      real(RK), intent(inout) :: dest(d, n)
      integer(IK)             :: i, j, fixed_indices(n-f)
!
      do concurrent(j=1:f, i=1:d)
        dest(i, j) = source(i, free_indices(j))
      end do
!
      call calc_fixed_indices(n, f, free_indices, fixed_indices)
!
      do concurrent(j=1:n - f, i=1:d)
        dest(i, j + f) = source(i, fixed_indices(j))
      end do
!
    end subroutine pack_matrix
!
  end subroutine lower_bound
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
    integer(IK)                       :: f
!
    if (d < 1 .or. n < 1) then
      res = 0
    else
      f = SIZE(free_m_indices)
      res = 2 + 3 * d * n + MAX(d * d * 2 + Kabsch_worksize(d), f * f * 2 + procrustes_worksize(f))
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
    integer(IK)                       :: i, j, maxiter_
    real(RK)                          :: threshold_
    integer(IK)                       :: f, g, mn, dd, df, ff, dfg, dmn
    integer(IK)                       :: conv, prev
    integer(IK)                       :: ix, iy, iz, drot, nrot, bs
    integer(IK)                       :: iw1, iw2
!
    f = SIZE(free_m_indices)
    g = SIZE(free_n_indices)
!
    if (f < 1 .or. g < 1) return
!
    mn = m * n
    dd = d * d
    df = d * f
    ff = f * f
    dmn = d * m * n
    dfg = d * f * g
    bs = ff + Procrustes_worksize(f)
!
    conv = 1
    prev = conv + 1
    ix   = prev + 1   ! copy of X
    iy   = ix   + dmn ! copy of Y
    iz   = iy   + dmn ! Z = P@Y
    drot = iz   + dmn ! rotation matrix P
    iw1  = drot + dd
    nrot = iz   + dmn ! rotation matrix Q
    iw2  = nrot + ff * g
!
    maxiter_   = optarg(maxiter,   DEF_maxiter)
    threshold_ = optarg(threshold, DEF_threshold) ** 2 * n
!
    call pack_matrix(d, m, n, f, g, free_m_indices, free_n_indices, X, w(ix))
    call pack_matrix(d, m, n, f, g, free_m_indices, free_n_indices, Y, w(iy))
!
    do concurrent(j=0:dmn-1)
      w(iz + j) = w(iy + j)
    end do
!
    w(conv) = sd(d, mn, w(ix), w(iy))
!
    do i = 1, maxiter_
!
      w(prev) = w(conv)
!
      call R_rotation(d, mn, w(ix), w(iy), w(iz), w(drot), w(iw1))
!print'(3f9.3)', w(iz:iz+dmn-1)
!print*
!
      do concurrent(j=0:g-1)
        call block_P_rotation(d, f, w(ix + j * df), w(iz + j * df), w(nrot + j * ff), w(iw2 + j * bs))
      end do
!
!print'(3f9.3)', w(iz:iz+dmn-1)
!print*
      w(conv) = sd(d, mn, w(ix), w(iz))
      if (ABS(w(prev) - w(conv)) < threshold_) exit
!
      do concurrent(j=0:g-1)
        call DGEMM('N', 'T', d, f, f, ONE, w(iy + j * df), d, w(nrot + j * ff), f, ZERO, w(iz + j * df), d)
      enddo

      do concurrent(j=dfg:dmn-1)
        w(iz + j) = w(iy + j)
      enddo
!
    enddo
!
    w(conv) = SQRT(w(conv) / n)
!
  contains
!
    pure subroutine pack_matrix(d, m, n, f, g, free_m_indices, free_n_indices, source, dest)
      integer(IK), intent(in) :: d, m, n, f, g, free_m_indices(f), free_n_indices(g)
      real(RK), intent(in)    :: source(*)
      real(RK), intent(inout) :: dest(*)
      integer(IK)             :: ifx
!
      call pack_free_matrix(d, m, n, f, g, free_m_indices, free_n_indices, source, dest)
      ifx = d * f * g + 1
      call pack_fixed_matrix(d, m, n, f, g, free_m_indices, free_n_indices, source, dest(ifx))
!
    end subroutine pack_matrix
!
    pure subroutine pack_free_matrix(d, m, n, f, g, free_m_indices, free_n_indices, source, dest)
      integer(IK), intent(in) :: d, m, n, f, g, free_m_indices(f), free_n_indices(g)
      real(RK), intent(in)    :: source(d, m, n)
      real(RK), intent(inout) :: dest(d, f, g)
      integer(IK)             :: i, j, k
!
      do concurrent(k=1:g, j=1:f, i=1:d)
        dest(i, j, k) = source(i, free_m_indices(j), free_n_indices(k))
      end do
!
    end subroutine pack_free_matrix
!
    pure subroutine pack_fixed_matrix(d, m, n, f, g, free_m_indices, free_n_indices, source, dest)
      integer(IK), intent(in) :: d, m, n, f, g, free_m_indices(f), free_n_indices(g)
      real(RK), intent(in)    :: source(d, m, n)
      real(RK), intent(inout) :: dest(*)
      integer(IK)             :: fixed_m_indices(m - f), fixed_n_indices(n - g)
      integer(IK)             :: m_f, n_g, dm, dm_f, ip
      integer(IK)             :: i, j, k
!
      call calc_fixed_indices(m, f, free_m_indices, fixed_m_indices)
      call calc_fixed_indices(n, g, free_n_indices, fixed_n_indices)
!
      m_f = m - f
      dm_f = d * m_f
      do concurrent(k=1:g, j=1:m_f, i=1:d)
        block
          integer(IK) :: p
          p = i + (j - 1) * d + (k - 1) * dm_f
          dest(p) = source(i, fixed_m_indices(j), free_n_indices(k))
        end block
      end do
!
      ip = g * dm_f
      n_g = n - g
      dm = d * m
!
      do concurrent(k=1:n_g, j=1:m, i=1:d)
        block
          integer(IK) :: p
          p = ip + i + (j - 1) * d + (k - 1) * dm
          dest(p) = source(i, j, fixed_n_indices(k))
        end block
      end do
!
    end subroutine pack_fixed_matrix
!
  end subroutine block_lower_bound
!
  pure subroutine R_rotation(d, n, x, y, z, R, w)
    integer(IK), intent(in)           :: d, n
    real(RK), intent(in)              :: x(*), y(*)
    real(RK), intent(inout)           :: z(*), R(*), w(*)
    integer(IK)                       :: iw
!
      iw = d * d + 1
!
!!    M = X@Z^T
      call DGEMM('N', 'T', d, d, n, ONE, x, d, z, d, ZERO, w, d)
!!    calc R from M
      call Kabsch(d, w, R, w(iw))
!!    Z = R@Y
      call DGEMM('N', 'N', d, n, d, ONE, R, d, y, d, ZERO, z, d)
!
  end subroutine R_rotation
!
  pure subroutine block_P_rotation(d, f, x, z, P, w)
    integer(IK), intent(in)           :: d, f
    real(RK), intent(in)              :: x(*)
    real(RK), intent(inout)           :: z(*), P(*), w(*)
    integer(IK)                       :: i, iw, df
!
    df = d * f
    iw = f * f + 1
!
!!  M = X^T@Z
    call DGEMM('T', 'N', f, f, d, ONE, x, d, z, d, ZERO, w, f)
!!  calc P from M
    call Procrustes(f, w, P, w(iw))
!!  Z' = Z@P
    call DGEMM('N', 'T', d, f, f, ONE, z, d, P, f, ZERO, w(iw), d)
!!  Z = Z'
    iw = iw - 1
    do concurrent(i=1:df)
      z(i) = w(iw + i)
    end do
!
  end subroutine block_P_rotation
!
  pure subroutine calc_fixed_indices(n, f, free_indices, res)
    integer(IK), intent(in)    :: n, f, free_indices(f)
    integer(IK), intent(inout) :: res(n-f)
    integer(IK)                :: i, j
    j = 1
    do i = 1, n
      if (ANY(i == free_indices)) cycle
      res(j) = i
      if (j == f) return
      j = j + 1
    end do
  end subroutine calc_fixed_indices
!
end module mod_lower_bound
