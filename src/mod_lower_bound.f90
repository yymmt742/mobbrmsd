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
  public :: lower_bound_worksize, lower_bound
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
    integer(IK)                       :: conv, prev
    integer(IK)                       :: ix, iy, iz, dcov, drot, ncov, nrot
    integer(IK)                       :: iw1, iw2, iw3
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
    dcov = drot + dd ! covariance matrix, XY^T
    nrot = iz   + dn ! rotation matrix Q
    ncov = nrot + ff ! covariance matrix, X^TY
    iw1  = dcov + dd
    iw2  = nrot + ff
    iw3  = iw2  + df
!
    maxiter_   = optarg(maxiter,   1000_IK)
    threshold_ = optarg(threshold, 1.0E-12_RK) ** 2 * n
!
    do concurrent(j=1:dn)
      w(ix + j - 1) = X(j)
    end do
    do concurrent(j=1:dn)
      w(iy + j - 1) = Y(j)
    end do
    do concurrent(j=1:dn)
      w(iz + j - 1) = w(iy + j - 1)
    end do
!
    w(conv) = sd(d, n, w(ix), w(iy))
!
    do i = 1, maxiter_
!
      w(prev) = w(conv)
!
      call cov(d, n, w(ix:ix+dn-1), w(iz:iz+dn-1), w(dcov:dcov+dd-1))
      call Kabsch(d, n, w(dcov), w(drot), w(iw1))
      call R_rotation(d, n, w(iy), w(drot), w(iz))
!
      call cov_row_major(d, free_indices, w(ix:ix+dn-1), w(iz:iz+dn-1), w(ncov:ncov+ff-1))
      call Procrustes(f, w(ncov), w(nrot), w(iw2))
      call partial_P_rotation(d, n, f, free_indices, w(nrot), w(iz), w(iw2), w(iw3))
!
      w(conv) = sd(d, n, w(ix), w(iz))
      if (ABS(w(prev) - w(conv)) < threshold_) exit
!
      do concurrent(j=1:dn)
        w(iz + j - 1) = w(iy + j - 1)
      enddo
!
      call partial_P_rotation(d, n, f, free_indices, w(nrot), w(iz), w(iw2), w(iw3))
!
    enddo
!
    w(conv) = SQRT(w(conv) / n)
!
  end subroutine lower_bound
!
  pure subroutine R_rotation(d, n, Y, R, Z)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: Y(d, n), R(d, d)
    real(RK), intent(inout) :: Z(d, n)
    Z = MATMUL(R, Y)
  end subroutine R_rotation

  pure subroutine partial_P_rotation(d, n, f, free_indices, P, Z, W1, W2)
    integer(IK), intent(in) :: d, n, f, free_indices(:)
    real(RK), intent(in)    :: P(f, f)
    real(RK), intent(inout) :: Z(d, n), W1(d, f), W2(d, f)
    integer(IK)             :: i
!
    do concurrent(i=1:f)
      W1(:, i) = Z(:, free_indices(i))
    enddo
!
    W2 = MATMUL(W1, TRANSPOSE(P))
!
    do concurrent(i=1:f)
      Z(:, free_indices(i)) = W2(:, i)
    enddo
!
  end subroutine partial_P_rotation
!
end module mod_lower_bound
