!> \brief \b DLASQ5 computes one dqds transform in ping-pong form. Used by sbdsqr and sstegr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASQ5 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq5.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq5.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq5.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN,
!                          DNM1, DNM2, IEEE, EPS )
!
!       .. Scalar Arguments ..
!       LOGICAL            IEEE
!       INTEGER            I0, N0, PP
!       DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU, SIGMA, EPS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASQ5 computes one dqds transform in ping-pong form, one
!> version for IEEE machines another for non IEEE machines.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] I0
!> \verbatim
!>          I0 is INTEGER
!>        First index.
!> \endverbatim
!>
!> \param[in] N0
!> \verbatim
!>          N0 is INTEGER
!>        Last index.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension ( 4*N )
!>        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
!>        an extra argument.
!> \endverbatim
!>
!> \param[in] PP
!> \verbatim
!>          PP is INTEGER
!>        PP=0 for ping, PP=1 for pong.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>        This is the shift.
!> \endverbatim
!>
!> \param[in] SIGMA
!> \verbatim
!>          SIGMA is DOUBLE PRECISION
!>        This is the accumulated shift up to this step.
!> \endverbatim
!>
!> \param[out] DMIN
!> \verbatim
!>          DMIN is DOUBLE PRECISION
!>        Minimum value of d.
!> \endverbatim
!>
!> \param[out] DMIN1
!> \verbatim
!>          DMIN1 is DOUBLE PRECISION
!>        Minimum value of d, excluding D( N0 ).
!> \endverbatim
!>
!> \param[out] DMIN2
!> \verbatim
!>          DMIN2 is DOUBLE PRECISION
!>        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!> \endverbatim
!>
!> \param[out] DN
!> \verbatim
!>          DN is DOUBLE PRECISION
!>        d(N0), the last value of d.
!> \endverbatim
!>
!> \param[out] DNM1
!> \verbatim
!>          DNM1 is DOUBLE PRECISION
!>        d(N0-1).
!> \endverbatim
!>
!> \param[out] DNM2
!> \verbatim
!>          DNM2 is DOUBLE PRECISION
!>        d(N0-2).
!> \endverbatim
!>
!> \param[in] IEEE
!> \verbatim
!>          IEEE is LOGICAL
!>        Flag for IEEE or non IEEE arithmetic.
!> \endverbatim
!>
!> \param[in] EPS
!> \verbatim
!>          EPS is DOUBLE PRECISION
!>        This is the value of epsilon used.
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
pure subroutine DLASQ5(I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, &
               &       DN, DNM1, DNM2, IEEE, EPS)
! use LA_CONSTANTS, only: RK => dp
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  logical, intent(in)   :: IEEE
  integer, intent(in)   :: I0, N0, PP
  real(RK), intent(in)  :: SIGMA, EPS
  real(RK), intent(out) :: TAU, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: Z(*)
!     ..
!  =====================================================================
!     .. Parameter ..
! real(RK), parameter     :: ZERO = 0.0_RK
! real(RK), parameter     :: HALF = 0.5_RK
!     ..
!     .. Local Scalars ..
  integer                :: J4, J4P2
  real(RK)               :: D, EMIN, TEMP, DTHRESH
!     ..
!     .. Intrinsic Functions ..
  intrinsic              :: MIN
!     ..
!     .. Executable Statements ..
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
!        Code for IEEE arithmetic.
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
!        Unroll last two steps.
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
!        Code for non IEEE arithmetic.
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
!        Unroll last two steps.
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
  else
!     This is the version that sets d's to zero if they are small enough
    J4 = 4 * I0 + PP - 3
    EMIN = Z(J4 + 4)
    D = Z(J4) - TAU
    DMIN = D
    DMIN1 = -Z(J4)
    if (IEEE) then
!
!     Code for IEEE arithmetic.
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
!     Unroll last two steps.
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
!     Code for non IEEE arithmetic.
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
!     Unroll last two steps.
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
  end if
!
  Z(J4 + 2) = DN
  Z(4 * N0 - PP) = EMIN
  return
!
!     End of DLASQ5
!
end subroutine DLASQ5

