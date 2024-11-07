!| mobbrmsd_SLASQ3 checks for deflation, computes a shift and calls dqds. Used by sbdsqr.
!
!  mobbrmsd_SLASQ3 checks for deflation, computes a shift (TAU) and calls dqds.
!  In case of failure it changes shifts, and tries again until output
!  is positive.
!
!  These are passed as arguments in order to save their values
!  between calls to mobbrmsd_SLASQ3.
!
!  Reference SLASQ3 is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_SLASQ3(I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, &
               &                NFAIL, ITER, NDIV, IEEE, TTYPE, DMIN1, &
               &                DMIN2, DN, DN1, DN2, G, TAU)
  implicit none
  integer, intent(in)     :: I0
!!  First index.
!!
  integer, intent(inout)  :: N0
!!  Last index.
!!
  real(RK), intent(inout) :: Z(*)
!!  REAL array, dimension ( 4*N0 )
!!  Z holds the qd array.
!!
  integer, intent(inout)  :: PP
!!  PP=0 for ping, PP=1 for pong.
!!  PP=2 indicates that flipping was applied to the Z array
!!  and that the initial tests for deflation should not be
!!  performed.
!!
  real(RK), intent(out)   :: DMIN
!!  Minimum value of d.
!!
  real(RK), intent(out)   :: SIGMA
!!  Sum of shifts used in current segment.
!!
  real(RK), intent(inout) :: DESIG
!!  Lower order part of SIGMA
!!
  real(RK), intent(inout) :: QMAX
!!  Maximum value of q.
!!
  integer, intent(inout)  :: NFAIL
!!  Increment NFAIL by 1 each time the shift was too big.
!!
  integer, intent(inout)  :: ITER
!!  Increment ITER by 1 for each iteration.
!!
  integer, intent(inout)  :: NDIV
!!  Increment NDIV by 1 for each division.
!!
  logical, intent(in)     :: IEEE
!!   IEEE is LOGICAL
!!  Flag for IEEE or non IEEE arithmetic (passed to mobbrmsd_DLASQ5).
!!
  integer, intent(inout)  :: TTYPE
!!  Shift type.
!!
  real(RK), intent(inout) :: DMIN1
!!  REAL
!!
  real(RK), intent(inout) :: DMIN2
!!  REAL
!!
  real(RK), intent(inout) :: DN
!!  REAL
!!
  real(RK), intent(inout) :: DN1
!!  REAL
!!
  real(RK), intent(inout) :: DN2
!!  REAL
!!
  real(RK), intent(inout) :: G
!!  REAL
!!
  real(RK), intent(inout) :: TAU
!!  REAL
!!
  real(RK), parameter :: CBIAS = 1.50E0
  integer :: IPN4, J4, N0IN, NN
  real    :: EPS, S, T, TEMP, TOL, TOL2
  intrinsic :: ABS, MAX, MIN, SQRT
! interface
!   include 'sisnan.h'
!   include 'slamch.h'
!   include 'slasq4.h'
!   include 'slasq5.h'
!   include 'slasq6.h'
! end interface
!
  N0IN = N0
  EPS = mobbrmsd_SLAMCH('Precision')
  TOL = EPS * HUNDRD
  TOL2 = TOL**2
!
! Check for deflation.
!
10 continue
!
  if (N0 < I0) return
  if (N0 == I0) GO TO 20
  NN = 4 * N0 + PP
  if (N0 == (I0 + 1)) GO TO 40
!
! Check whether E(N0 - 1) is negligible, 1 eigenvalue.
!
  if (Z(NN - 5) > TOL2 * (SIGMA + Z(NN - 3)) .and. Z(NN - 2 * PP - 4) > TOL2 * Z(NN - 7)) GO TO 30
!
20 continue
!
  Z(4 * N0 - 3) = Z(4 * N0 + PP - 3) + SIGMA
  N0 = N0 - 1
  GO TO 10
!
! Check whether E(N0 - 2) is negligible, 2 eigenvalues.
!
30 continue
!
  if (Z(NN - 9) > TOL2 * SIGMA .and. Z(NN - 2 * PP - 8) > TOL2 * Z(NN - 11)) GO TO 50
!
40 continue
!
  if (Z(NN - 3) > Z(NN - 7)) then
    S = Z(NN - 3)
    Z(NN - 3) = Z(NN - 7)
    Z(NN - 7) = S
  end if
  T = HALF * ((Z(NN - 7) - Z(NN - 3)) + Z(NN - 5))
  if (Z(NN - 5) > Z(NN - 3) * TOL2 .and. T /= ZERO) then
    S = Z(NN - 3) * (Z(NN - 5) / T)
    if (S <= T) then
      S = Z(NN - 3) * (Z(NN - 5) / (T * (ONE + SQRT(ONE + S / T))))
    else
      S = Z(NN - 3) * (Z(NN - 5) / (T + SQRT(T) * SQRT(T + S)))
    end if
    T = Z(NN - 7) + (S + Z(NN - 5))
    Z(NN - 3) = Z(NN - 3) * (Z(NN - 7) / T)
    Z(NN - 7) = T
  end if
  Z(4 * N0 - 7) = Z(NN - 7) + SIGMA
  Z(4 * N0 - 3) = Z(NN - 3) + SIGMA
  N0 = N0 - 2
  GO TO 10
!
50 continue
  if (PP == 2) PP = 0
!
!Reverse the qd - array, if warranted.
!
  if (DMIN <= ZERO .or. N0 < N0IN) then
    if (CBIAS * Z(4 * I0 + PP - 3) < Z(4 * N0 + PP - 3)) then
      IPN4 = 4 * (I0 + N0)
      do J4 = 4 * I0, 2 * (I0 + N0 - 1), 4
        TEMP = Z(J4 - 3)
        Z(J4 - 3) = Z(IPN4 - J4 - 3)
        Z(IPN4 - J4 - 3) = TEMP
        TEMP = Z(J4 - 2)
        Z(J4 - 2) = Z(IPN4 - J4 - 2)
        Z(IPN4 - J4 - 2) = TEMP
        TEMP = Z(J4 - 1)
        Z(J4 - 1) = Z(IPN4 - J4 - 5)
        Z(IPN4 - J4 - 5) = TEMP
        TEMP = Z(J4)
        Z(J4) = Z(IPN4 - J4 - 4)
        Z(IPN4 - J4 - 4) = TEMP
      end do
      if (N0 - I0 <= 4) then
        Z(4 * N0 + PP - 1) = Z(4 * I0 + PP - 1)
        Z(4 * N0 - PP) = Z(4 * I0 - PP)
      end if
      DMIN2 = MIN(DMIN2, Z(4 * N0 + PP - 1))
      Z(4 * N0 + PP - 1) = MIN(Z(4 * N0 + PP - 1), Z(4 * I0 + PP - 1), Z(4 * I0 + PP + 3))
      Z(4 * N0 - PP) = MIN(Z(4 * N0 - PP), Z(4 * I0 - PP), Z(4 * I0 - PP + 4))
      QMAX = MAX(QMAX, Z(4 * I0 + PP - 3), Z(4 * I0 + PP + 1))
      DMIN = -ZERO
    end if
  end if
!
! Choose a shift.
!
  call mobbrmsd_SLASQ4(I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, DN2, TAU, TTYPE, G)
!
! call dqds until DMIN > 0.
!
70 continue
!
  call mobbrmsd_SLASQ5(I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, DN1, DN2, IEEE, EPS)
!
  NDIV = NDIV + (N0 - I0 + 2)
  ITER = ITER + 1
!
! Check status.
!
  if (DMIN >= ZERO .and. DMIN1 >= ZERO) then
!
! Success.
!
    GO TO 90
!
  else if (DMIN < ZERO .and. DMIN1 > ZERO .and. &
        & Z(4 * (N0 - 1) - PP) < TOL * (SIGMA + DN1) .and. &
        & ABS(DN) < TOL * SIGMA) then
!
! Convergence hidden by negative DN.
!
    Z(4 * (N0 - 1) - PP + 2) = ZERO
    DMIN = ZERO
    GO TO 90
  else if (DMIN < ZERO) then
!
! TAU too big.select new TAU and try again.
!
    NFAIL = NFAIL + 1
    if (TTYPE < -22) then
!
! Failed twice. Play it safe.
!
      TAU = ZERO
    else if (DMIN1 > ZERO) then
!
! Late failure.Gives excellent shift.
!
      TAU = (TAU + DMIN) * (ONE - TWO * EPS)
      TTYPE = TTYPE - 11
    else
!
! Early failure.Divide by 4.
!
      TAU = QURTR * TAU
      TTYPE = TTYPE - 12
    end if
    GO TO 70
  else if (IEEE_IS_NAN(DMIN)) then
!
! NaN.
!
    if (TAU == ZERO) then
      GO TO 80
    else
      TAU = ZERO
      GO TO 70
    end if
  else
!
! Possible underflow.Play it safe.
!
    GO TO 80
  end if
!
! Risk of underflow.
!
80 continue
  call mobbrmsd_SLASQ6(I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2)
  NDIV = NDIV + (N0 - I0 + 2)
  ITER = ITER + 1
  TAU = ZERO
!
90 continue
  if (TAU < SIGMA) then
    DESIG = DESIG + TAU
    T = SIGMA + DESIG
    DESIG = DESIG - (T - SIGMA)
  else
    T = SIGMA + TAU
    DESIG = SIGMA - (T - TAU) + DESIG
  end if
  SIGMA = T
!
  return
!
!end of mobbrmsd_SLASQ3
!
end

