!| Calculate min_{R,P} |X-R@Y@P|^2 for R\in\R^{d,d} and P\in \R^{n,n}.
!  Here, R@R^T=I, det(R)=1, P@P^T=I are satisfied.
module mod_lower_bound
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_optarg
  use mod_rmsd
  use mod_Kabsch
  use mod_Procrustes
  implicit none
  private
  public :: lower_bound_worksize, lower_bound
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
    integer(IK)                       :: n_indices(n)
    integer(IK)                       :: conv, prev
    integer(IK)                       :: ix, iy, iz, rcov, rrot, ncov, nrot
    integer(IK)                       :: iw1, iw2
!
    f = SIZE(free_indices)
    if (f < 1) then
      w(1) = rmsd(d, n, X, Y)
      return
    end if
!
    dn = d * n
    dd = d * d
    df = d * f
    ff = f * f
!&<
    block
      integer(IK) :: l, u
      u = 0
      l = u + 1 ; u = u + 1  ; conv = l
      l = u + 1 ; u = u + 1  ; prev = l
      l = u + 1 ; u = u + dn ; ix   = l ! copy of X
      l = u + 1 ; u = u + dn ; iy   = l ! copy of Y
      l = u + 1 ; u = u + dn ; iz   = l ! Z = P@Y
      l = u + 1 ; u = u + dd ; rrot = l ! rotation matrix R
      l = u + 1 ; u = u + dd ; rcov = l ! covariance matrix, X^TY
      l = u + 1 ; u = rrot   ; iw1  = l
      l = u + 1 ; u = u + ff ; nrot = l ! rotation matrix Q
      l = u + 1 ; u = u + ff ; ncov = l ! covariance matrix, X^TY
      l = u + 1 ; u = u + 1  ; iw2  = l
    end block
!>&
    maxiter_   = optarg(maxiter,   DEF_maxiter)
    threshold_ = optarg(threshold, DEF_threshold) ** 2 * n / 2
!
    call calc_swap_indices(n, f, free_indices, n_indices)
!
    call pack_matrix(d, n, n_indices, X, w(ix))
    call pack_matrix(d, n, n_indices, Y, w(iy))
!
    do concurrent(j=0:dn - 1)
      w(iz + j) = w(iy + j)
    end do
!
    w(conv) = -RHUGE
!
   do i = 1, maxiter_
!
      w(prev) = w(conv)
!
!!!!  P rotation  !!!!
!
!!    M = X^T@Z
      call DGEMM('T', 'N', f, f, d, ONE, w(ix), d, w(iz), d, ZERO, w(ncov), f)
!!    calc P from M
      call Procrustes(f, w(ncov), w(nrot), w(iw2))
!
!!    Z = Y@P^T
      call DGEMM('N', 'T', d, f, f, ONE, w(iy), d, w(nrot), f, ZERO, w(iz), d)
      do concurrent(j=df:dn - 1)
        w(iz + j) = w(iy + j)
      end do
!
!!!!  R rotation  !!!!
!
!!    M = X@Z^T
      call DGEMM('N', 'T', d, d, n, ONE, w(ix), d, w(iz), d, ZERO, w(rcov), d)
!!    get optimal R
      call Kabsch(d, w(rcov), w(rrot), w(iw1))
!
!!    trace(X@R@Y@P) = trace(R^T@(Y@P@X))
!
      w(conv) = ZERO
      do j = 0, dd - 1
        w(conv) = w(conv) + w(rcov + j) * w(rrot + j)
      end do
!
      if (ABS(w(prev) - w(conv)) < threshold_) then
!!    Y = R@Z
        call DGEMM('N', 'N', d, n, d, ONE, w(rrot), d, w(iz), d, ZERO, w(iy), d)
        exit
      else
!!    Z = R@Y
        call DGEMM('N', 'N', d, n, d, ONE, w(rrot), d, w(iy), d, ZERO, w(iz), d)
      end if
!
    enddo
!
    w(conv) = rmsd(d, n, w(ix), w(iy))
!
  contains
!
    pure subroutine pack_matrix(d, n, n_indices, source, dest)
      integer(IK), intent(in) :: d, n, n_indices(n)
      real(RK), intent(in)    :: source(d, n)
      real(RK), intent(inout) :: dest(d, n)
      integer(IK)             :: i, j
!
      do concurrent(j=1:n, i=1:d)
        dest(i, j) = source(i, n_indices(j))
      end do
!
    end subroutine pack_matrix
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
  end subroutine lower_bound
!
end module mod_lower_bound
