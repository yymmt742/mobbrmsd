!| mobbrmsd_SLASQ5 computes one dqds transform in ping-pong form.
!
!  mobbrmsd_SLASQ5 computes one dqds transform in ping-pong form,
!  one version for IEEE machines another for non IEEE machines.
!
!  Reference SLASQ5 is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
pure subroutine mobbrmsd_SLASQ5(I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, &
               &                DMIN2, DN, DNM1, DNM2, IEEE, EPS)
  implicit none
  integer, intent(in)   :: I0
!!  First index.
!!
  integer, intent(in)   :: N0
!!  Last index.
!!
  real(RK), intent(inout) :: Z(*)
!!  REAL array, dimension ( 4*N )
!!  Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
!!  an extra argument.
!!
  integer, intent(in)   :: PP
!!  PP=0 for ping, PP=1 for pong.
!!
  real(RK), intent(out) :: TAU
!!  This is the shift.
!!
  real(RK), intent(in)  :: SIGMA
!!  This is the accumulated shift up to this step.
!!
  real(RK), intent(out) :: DMIN
!!  Minimum value of d.
!!
  real(RK), intent(out) :: DMIN1
!!  Minimum value of d, excluding D( N0 ).
!!
  real(RK), intent(out) :: DMIN2
!!  Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!!
  real(RK), intent(out) :: DN
!!  d(N0), the last value of d.
!!
  real(RK), intent(out) :: DNM1
!!  d(N0-1).
!!
  real(RK), intent(out) :: DNM2
!!  d(N0-2).
!!
  logical, intent(in)   :: IEEE
!!  Flag for IEEE or non IEEE arithmetic.
!!
  real(RK), intent(in)  :: EPS
!!  This is the value of epsilon used.
!!
  integer :: J4, J4P2
  real(RK) :: D, EMIN, TEMP, DTHRESH
  intrinsic :: MIN
!
  if ((N0 - I0 - 1) <= 0) return
!
  DTHRESH = EPS * (SIGMA + TAU)
  if (TAU < DTHRESH * HALF) TAU = ZERO
  if (TAU /= ZERO) then
    J4 = 4 * I0 + PP - 3
    EMIN = Z(J4 + 4)
    D = Z(J4) - TAU
    DMIN = D
    DMIN1 = -Z(J4)
!
    if (IEEE) then
!
! Code for IEEE arithmetic.
!
      if (PP == 0) then
        do J4 = 4 * I0, 4 * (N0 - 3), 4
          Z(J4 - 2) = D + Z(J4 - 1)
          TEMP = Z(J4 + 1) / Z(J4 - 2)
          D = D * TEMP - TAU
          DMIN = MIN(DMIN, D)
          Z(J4) = Z(J4 - 1) * TEMP
          EMIN = MIN(Z(J4), EMIN)
        end do
      else
        do J4 = 4 * I0, 4 * (N0 - 3), 4
          Z(J4 - 3) = D + Z(J4)
          TEMP = Z(J4 + 2) / Z(J4 - 3)
          D = D * TEMP - TAU
          DMIN = MIN(DMIN, D)
          Z(J4 - 1) = Z(J4) * TEMP
          EMIN = MIN(Z(J4 - 1), EMIN)
        end do
      end if
!
! Unroll last two steps.
!
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4 * (N0 - 2) - PP
      J4P2 = J4 + 2 * PP - 1
      Z(J4 - 2) = DNM2 + Z(J4P2)
      Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
      DNM1 = Z(J4P2 + 2) * (DNM2 / Z(J4 - 2)) - TAU
      DMIN = MIN(DMIN, DNM1)
!
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2 * PP - 1
      Z(J4 - 2) = DNM1 + Z(J4P2)
      Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
      DN = Z(J4P2 + 2) * (DNM1 / Z(J4 - 2)) - TAU
      DMIN = MIN(DMIN, DN)
!
    else
!
! Code for non IEEE arithmetic.
!
      if (PP == 0) then
        do J4 = 4 * I0, 4 * (N0 - 3), 4
          Z(J4 - 2) = D + Z(J4 - 1)
          if (D < ZERO) then
            return
          else
            Z(J4) = Z(J4 + 1) * (Z(J4 - 1) / Z(J4 - 2))
            D = Z(J4 + 1) * (D / Z(J4 - 2)) - TAU
          end if
          DMIN = MIN(DMIN, D)
          EMIN = MIN(EMIN, Z(J4))
        end do
      else
        do J4 = 4 * I0, 4 * (N0 - 3), 4
          Z(J4 - 3) = D + Z(J4)
          if (D < ZERO) then
            return
          else
            Z(J4 - 1) = Z(J4 + 2) * (Z(J4) / Z(J4 - 3))
            D = Z(J4 + 2) * (D / Z(J4 - 3)) - TAU
          end if
          DMIN = MIN(DMIN, D)
          EMIN = MIN(EMIN, Z(J4 - 1))
        end do
      end if
!
! Unroll last two steps.
!
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4 * (N0 - 2) - PP
      J4P2 = J4 + 2 * PP - 1
      Z(J4 - 2) = DNM2 + Z(J4P2)
      if (DNM2 < ZERO) then
        return
      else
        Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
        DNM1 = Z(J4P2 + 2) * (DNM2 / Z(J4 - 2)) - TAU
      end if
      DMIN = MIN(DMIN, DNM1)
!
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2 * PP - 1
      Z(J4 - 2) = DNM1 + Z(J4P2)
      if (DNM1 < ZERO) then
        return
      else
        Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
        DN = Z(J4P2 + 2) * (DNM1 / Z(J4 - 2)) - TAU
      end if
      DMIN = MIN(DMIN, DN)
!
    end if
!
  else
! This is the version that sets d's to zero if they are small enough
    J4 = 4 * I0 + PP - 3
    EMIN = Z(J4 + 4)
    D = Z(J4) - TAU
    DMIN = D
    DMIN1 = -Z(J4)
    if (IEEE) then
!
! Code for IEEE arithmetic.
!
      if (PP == 0) then
        do J4 = 4 * I0, 4 * (N0 - 3), 4
          Z(J4 - 2) = D + Z(J4 - 1)
          TEMP = Z(J4 + 1) / Z(J4 - 2)
          D = D * TEMP - TAU
          if (D < DTHRESH) D = ZERO
          DMIN = MIN(DMIN, D)
          Z(J4) = Z(J4 - 1) * TEMP
          EMIN = MIN(Z(J4), EMIN)
        end do
      else
        do J4 = 4 * I0, 4 * (N0 - 3), 4
          Z(J4 - 3) = D + Z(J4)
          TEMP = Z(J4 + 2) / Z(J4 - 3)
          D = D * TEMP - TAU
          if (D < DTHRESH) D = ZERO
          DMIN = MIN(DMIN, D)
          Z(J4 - 1) = Z(J4) * TEMP
          EMIN = MIN(Z(J4 - 1), EMIN)
        end do
      end if
!
! Unroll last two steps.
!
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4 * (N0 - 2) - PP
      J4P2 = J4 + 2 * PP - 1
      Z(J4 - 2) = DNM2 + Z(J4P2)
      Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
      DNM1 = Z(J4P2 + 2) * (DNM2 / Z(J4 - 2)) - TAU
      DMIN = MIN(DMIN, DNM1)
!
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2 * PP - 1
      Z(J4 - 2) = DNM1 + Z(J4P2)
      Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
      DN = Z(J4P2 + 2) * (DNM1 / Z(J4 - 2)) - TAU
      DMIN = MIN(DMIN, DN)
!
    else
!
! Code for non IEEE arithmetic.
!
      if (PP == 0) then
        do J4 = 4 * I0, 4 * (N0 - 3), 4
          Z(J4 - 2) = D + Z(J4 - 1)
          if (D < ZERO) then
            return
          else
            Z(J4) = Z(J4 + 1) * (Z(J4 - 1) / Z(J4 - 2))
            D = Z(J4 + 1) * (D / Z(J4 - 2)) - TAU
          end if
          if (D < DTHRESH) D = ZERO
          DMIN = MIN(DMIN, D)
          EMIN = MIN(EMIN, Z(J4))
        end do
      else
        do J4 = 4 * I0, 4 * (N0 - 3), 4
          Z(J4 - 3) = D + Z(J4)
          if (D < ZERO) then
            return
          else
            Z(J4 - 1) = Z(J4 + 2) * (Z(J4) / Z(J4 - 3))
            D = Z(J4 + 2) * (D / Z(J4 - 3)) - TAU
          end if
          if (D < DTHRESH) D = ZERO
          DMIN = MIN(DMIN, D)
          EMIN = MIN(EMIN, Z(J4 - 1))
        end do
      end if
!
! Unroll last two steps.
!
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4 * (N0 - 2) - PP
      J4P2 = J4 + 2 * PP - 1
      Z(J4 - 2) = DNM2 + Z(J4P2)
      if (DNM2 < ZERO) then
        return
      else
        Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
        DNM1 = Z(J4P2 + 2) * (DNM2 / Z(J4 - 2)) - TAU
      end if
      DMIN = MIN(DMIN, DNM1)
!
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2 * PP - 1
      Z(J4 - 2) = DNM1 + Z(J4P2)
      if (DNM1 < ZERO) then
        return
      else
        Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
        DN = Z(J4P2 + 2) * (DNM1 / Z(J4 - 2)) - TAU
      end if
      DMIN = MIN(DMIN, DN)
!
    end if
!
  end if
  Z(J4 + 2) = DN
  Z(4 * N0 - PP) = EMIN
  return
!
! end of mobbrmsd_SLASQ5
!
end

