!| Calculate the rotation matrix that minimizes |X-RY|^2 using the Kabsch-Umeyama algorithm.
!  Here, RR^T=I and det(R)=1 are satisfied.
module mod_Kabsch
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_svd
  use mod_det
  implicit none
  private
  public :: Kabsch_worksize, Kabsch
  public :: get_rotation_matrix
!
  interface
    include 'dgesvd.h'
    include 'dsyev.h'
    include 'dgemm.h'
  end interface
!
contains
!
!| Calculate work array size for d*d matrix.
  pure elemental function Kabsch_worksize(d) result(res)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    integer(IK)                   :: res
!
    if (d < 1) then
      res = 0
    else
      res = svd_worksize(d) + d * d * 3 + d
    end if
!
  end function Kabsch_worksize
!
!| Calculate the rotation matrix from covariance matrix.
   subroutine Kabsch(d, cov, rot, w, ldcov)
    integer(IK), intent(in)       :: d
    !! matrix collumn dimension.
    real(RK), intent(in)          :: cov(*)
    !! target d*n array
    real(RK), intent(inout)       :: rot(*)
    !! rotation d*d matrix
    real(RK), intent(inout)       :: w(*)
    !! work array, must be larger than Kabsch_worksize(d)
    !! if row_major, must be larger than Kabsch_worksize(n)
    integer(IK), intent(in), optional :: ldcov
    integer(IK)                   :: dd, m, s, u, vt, iw
!
    if (d < 1) RETURN
    if (d == 1)then
      rot(1) = ONE
      RETURN
    endif
!
    dd = d * d
    m = 1
    u = m + dd
    vt = u + dd
    s = vt + dd
    iw = s + d
!
    w(:dd) = cov(:dd)
!
    call svd(d, w(m), w(s), w(u), w(vt), w(iw), ldx=ldcov)
!
    call DGEMM('N', 'N', d, d, d, ONE, w(u), d, w(vt), d, ZERO, w(s), d)
    call det_sign(d, w(s:s + dd - 1))
    if (w(s) < ZERO) w(u + dd - d:u + dd - 1) = -w(u + dd - d:u + dd - 1)
    call DGEMM('N', 'N', d, d, d, ONE, w(u), d, w(vt), d, ZERO, w(s), d)
!
    rot(:dd) = w(s:s+dd-1)
!
  end subroutine Kabsch
!
!| Calculate the rotation matrix from covariance matrix.
   subroutine get_rotation_matrix(d, n, X, Y, C, R, w, threshold)
    integer(IK), intent(in)        :: d, n
    !! matrix collumn dimension.
    real(RK), intent(in)           :: X(*)
    !! reference d*n array
    real(RK), intent(in)           :: Y(*)
    !! target d*n array
    real(RK), intent(in)           :: C(*)
    !! precomputed covariance matrix XY^T, d*d array
    real(RK), intent(inout)        :: R(*)
    !! rotation d*d matrix
    real(RK), intent(inout)        :: W(*)
    !! work array, must be larger than Kabsch_worksize(d)
    real(RK), intent(in), optional :: threshold
    !! threshold
    real(RK)                       :: threshold_
    integer(IK)                    :: z, s, u, cov, rot, iw, ds, dp
!
    if (d < 1) RETURN
    if (d == 1)then
      R(1) = ONE
      RETURN
    endif
!
    ds = MIN(d, n)
    u = 1
    s = u + d * d
    z = s + d
    iw = z + d * n
    call reduce_dimension(d, n, X, Y, w(z), w(u), w(s), w(iw))
!
print'(*(f9.3))',w(s:s+ds-1)
    threshold_ = 1.0E-8_RK
    if (PRESENT(threshold)) threshold_ = threshold
!
    dp = COUNT(w(s:s + ds - 1) > threshold_)
    if(dp<d)then
      call partial_Kabsch(d, n, X, Y, w(U), w(z), w(iw), R, w(iw+d*d))
      cov = s + ds
      iw = cov + dp * dp

!     call DGEMM('T', 'N', dp, d, d, ONE, w(u), d, C, d, ZERO, w(iw), d)
!     call DGEMM('N', 'N', dp, dp, d, ONE, w(iw), d, w(u), d, ZERO, w(cov), dp)
!     rot = cov + dp * dp
!     iw = rot + dp * dp
!print*,'R'
!print'(6f9.3)',R(:d*d)
!print*,'cov'
!print'(6f9.3)',w(cov:cov+dp*dp-1)
!      call Kabsch(dp, w(cov), w(rot), w(iw))
!print*,'rot'
!print'(6f9.3)',w(rot:rot+dp*dp-1)
!      call pack_rot(d, dp, n, w(rot), R)
!print*,'R'
!print'(6f9.3)',R(:d*d)
!print*,'U'
!print'(6f9.3)',w(u:u+d*d)
!      call DGEMM('T', 'N', d, d, d, ONE, w(u), d, R, d, ZERO, w(iw), d)
!      call DGEMM('N', 'N', d, d, d, ONE, w(iw), d, w(u), d, ZERO, R, d)
    else
      call Kabsch(d, C, R, w)
    endif
!
  end subroutine get_rotation_matrix
!
   subroutine reduce_dimension(d, n, X, Y, Z, U, S, W)
    integer(IK), intent(in)      :: d, n
    real(RK), intent(in)         :: X(d, n), Y(d, n)
    real(RK), intent(inout)      :: Z(d, n), U(d, d), S(d), W(*)
    integer(IK)                  :: i, j, lw, info
      do concurrent(i=1:d, j=1:n)
        Z(i, j) = X(i, j) - Y(i, j)
      end do
print*,'red,Z'
print'(6f9.3)', Z
print*
      call DGEMM('N', 'T', d, d, n,-ONE, Z, d, Z, d, ZERO, U, d)
      call DSYEV('V', 'L', d, U, d, S, w, -1, info)
      lw = NINT(w(1))
      call DSYEV('V', 'L', d, U, d, S, w, lw, info)
!     call DGESVD('A', 'N', d, n, Z, d, S, U, d, W, n, W, -1, info)
!     lw = NINT(w(1))
!     call DGESVD('A', 'N', d, n, Z, d, S, U, d, W, n, W(n*n+1), lw, info)
!U(:,:2) = -U(:,:2)
!U(:,:3) = 0D0
print*,'red,U'
print'(6f9.3)', U
print*
print*,'red,UU'
print'(6f9.3)', MATMUL(TRANSPOSE(U),U)
print*
S = -S
print*,'s'
print'(6f9.3)', S
print*,'U'
print'(6f9.3)', U
print*,'UZ'
print'(6f9.3)', SQRT(MATMUL(MATMUL(TRANSPOSE(U),X-Y),TRANSPOSE(MATMUL(TRANSPOSE(U),X-Y))))
print*

  end subroutine reduce_dimension
!
   subroutine partial_Kabsch(d, n, X, Y, U, Z, C, R, W)
    integer(IK), intent(in)      :: d, n
    real(RK), intent(in)         :: X(d, n), Y(d, n), U(d, d)
    real(RK), intent(inout)      :: Z(d, n), C(d, d), R(d, d), W(*)
    real(RK)                     :: WW(200)
    integer(IK)                  :: i, j, lw, info
    Z = MATMUL(U, X)
    C = MATMUL(MATMUL(Z, TRANSPOSE(Y)), TRANSPOSE(U))
    print*,'PC'
    print'(6f9.3)',C
!   call Kabsch(d, C, R, ww)
do concurrent(i=1:d,j=1:d)
R(i, j) = MERGE(1, 0, i == j)
enddo
    print*,'PR'
    print'(6f9.3)',R
    R = MATMUL(MATMUL(TRANSPOSE(U), R), U)
    print*,'PURU'
    print'(6f9.3)',R
    print*
  end subroutine partial_Kabsch
!
  pure subroutine pack_rot(d, dp, n, rot, R)
    integer(IK), intent(in) :: d, dp, n
    real(RK), intent(in)    :: rot(dp, dp)
    real(RK), intent(inout) :: R(d, d)
    integer(IK)             :: i, j
    do concurrent(i=1:d, j=1:d)
      if (i > dp .or. j > dp) then
        R(i, j) = MERGE(ONE, ZERO, i == j)
      else
        R(i, j) = rot(i, j)
      end if
    end do
  end subroutine pack_rot
!
end module mod_Kabsch
