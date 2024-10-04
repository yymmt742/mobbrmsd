!| Define spatial dimension, \(D\),
!  and provide optimized functions for dimension.
module mod_dimspec_functions
  use mod_kinds, only: RK, IK
#ifdef USE_REAL32
  use mod_mobbrmsd_lapack, only: SGEMM
#else
  use mod_mobbrmsd_lapack, only: DGEMM
#endif
  implicit none
  private
  public  :: D, DD, ND
  public  :: setup_dimension
  public  :: compute_com
  public  :: compute_cov
  public  :: covcopy
  !| Spatial dimension, \(D\).
  integer(IK), protected, save :: D = 3
  !| Square spatial dimension, \(D^2\).
  integer(IK), protected, save :: DD = 9
  !| Node memory size, defined by \(1 + 1 + D^2\).
  !  Let \([L, G, \mathbf{C}]\) be a node,
  !  where \(L, G\in\mathbb{R}\) and \(\mathbf{C}\in\mathbb{R}^{D\times D}\).
  integer(IK), protected, save :: ND = 9 + 2
!
  real(RK), parameter :: ZERO = 0.0_RK
  real(RK), parameter :: HALF = 0.5_RK
  real(RK), parameter :: ONETHIRD = 1.0_RK / 3.0_RK
  real(RK), parameter :: ONE = 1.0_RK
!
contains
!| Sets the dimensions of the space. <br>
!  Caution, this routine affects global.
  subroutine setup_dimension(d_)
    integer(IK), intent(in) :: d_
    D = MAX(1, d_)
    DD = D * D
    ND = DD + 2
  end subroutine setup_dimension
!
!| Calculate center of mass.
  pure subroutine compute_com(d, n, X, C)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *)
    real(RK), intent(inout) :: C(d)
    real(RK)                :: rn
    integer(IK)             :: i, j
    if (n < 0) then
      do concurrent(i=1:d)
        C(i) = ZERO
      end do
      return
    elseif (n == 1) then
      do concurrent(i=1:d)
        C(i) = X(i, 1)
      end do
      return
    elseif (n == 2) then
      do concurrent(i=1:d)
        C(i) = (X(i, 1) + X(i, 2)) * HALF
      end do
      return
    elseif (n == 3) then
      do concurrent(i=1:d)
        C(i) = (X(i, 1) + X(i, 2) + X(i, 3)) * ONETHIRD
      end do
      return
    end if
    do concurrent(i=1:d)
      C(i) = ZERO
    end do
    do j = 2, n, 2
      do concurrent(i=1:d)
        C(i) = C(i) + X(i, j - 1)
      end do
      do concurrent(i=1:d)
        C(i) = C(i) + X(i, j - 0)
      end do
    end do
    if (MODULO(n, 2) == 1) then
      do concurrent(i=1:d)
        C(i) = C(i) + X(i, n)
      end do
    end if
    rn = ONE / real(n, RK)
    do concurrent(i=1:d)
      C(i) = C(i) * rn
    end do
  end subroutine compute_com
!
!| Calculate covariance matrix.
  pure subroutine compute_cov(d, n, X, Y, C)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), Y(d, *)
    real(RK), intent(inout) :: C(d, d)
#ifdef USE_REAL32
    call SGEMM('N', 'T', D, D, n, ONE, Y, D, X, D, ZERO, C, D)
#else
    call DGEMM('N', 'T', D, D, n, ONE, Y, D, X, D, ZERO, C, D)
#endif
  end subroutine compute_cov

!| Compute \(\mathbf{Y} \gets \mathbf{X}-\bar{\mathbf{X}}\).
  pure subroutine covcopy(d, n, X, CX, Y)
    integer(IK), intent(in) :: d, n
    real(RK), intent(in)    :: X(d, *), CX(d)
    real(RK), intent(inout) :: Y(d, *)
    integer(IK)             :: i, j
    do concurrent(i=1:d, j=1:n)
      Y(i, j) = X(i, j) - CX(i)
    end do
  end subroutine covcopy
end module mod_dimspec_functions

