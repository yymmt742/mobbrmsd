!> \brief \b DLASQ6 computes one dqd transform in ping-pong form. Used by sbdsqr and sstegr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASQ6 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq6.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq6.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq6.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN,
!                          DNM1, DNM2 )
!
!       .. Scalar Arguments ..
!       INTEGER            I0, N0, PP
!       DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
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
!> DLASQ6 computes one dqd (shift equal to zero) transform in
!> ping-pong form, with protection against underflow and overflow.
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
!
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
pure subroutine DLASQ6(I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2)
! use LA_CONSTANTS, only: RK => dp
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)     :: I0, N0, PP
  real(RK), intent(out)   :: DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: Z(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
  integer  :: J4, J4P2
  real(RK) :: D, EMIN, SAFMIN, TEMP
!     ..
!     .. Intrinsic Functions ..
  intrinsic :: MIN
!
!     .. Parameter ..
! real(RK), parameter     :: ZERO = 0.0_RK
!     ..
!     .. External Function ..
! interface
!   pure elemental function DLAMCH(CMACH)
!     use LA_CONSTANTS, only: RK => dp
!     character(*), intent(in) :: CMACH
!     real(RK)                :: DLAMCH
!   end function DLAMCH
! end interface
!     ..
!     .. Executable Statements ..
!
  if ((N0 - I0 - 1) <= 0) return
!
  SAFMIN = DLAMCH('Safe minimum')
  J4 = 4 * I0 + PP - 3
  EMIN = Z(J4 + 4)
  D = Z(J4)
  DMIN = D
!
  if (PP == 0) then
    do J4 = 4 * I0, 4 * (N0 - 3), 4
      Z(J4 - 2) = D + Z(J4 - 1)
      if (Z(J4 - 2) == ZERO) then
        Z(J4) = ZERO
        D = Z(J4 + 1)
        DMIN = D
        EMIN = ZERO
      else if (SAFMIN * Z(J4 + 1) < Z(J4 - 2) .and. &
     &         SAFMIN * Z(J4 - 2) < Z(J4 + 1)) then
        TEMP = Z(J4 + 1) / Z(J4 - 2)
        Z(J4) = Z(J4 - 1) * TEMP
        D = D * TEMP
      else
        Z(J4) = Z(J4 + 1) * (Z(J4 - 1) / Z(J4 - 2))
        D = Z(J4 + 1) * (D / Z(J4 - 2))
      end if
      DMIN = MIN(DMIN, D)
      EMIN = MIN(EMIN, Z(J4))
    end do
  else
    do J4 = 4 * I0, 4 * (N0 - 3), 4
      Z(J4 - 3) = D + Z(J4)
      if (Z(J4 - 3) == ZERO) then
        Z(J4 - 1) = ZERO
        D = Z(J4 + 2)
        DMIN = D
        EMIN = ZERO
      else if (SAFMIN * Z(J4 + 2) < Z(J4 - 3) .and. &
     &         SAFMIN * Z(J4 - 3) < Z(J4 + 2)) then
        TEMP = Z(J4 + 2) / Z(J4 - 3)
        Z(J4 - 1) = Z(J4) * TEMP
        D = D * TEMP
      else
        Z(J4 - 1) = Z(J4 + 2) * (Z(J4) / Z(J4 - 3))
        D = Z(J4 + 2) * (D / Z(J4 - 3))
      end if
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
  if (Z(J4 - 2) == ZERO) then
    Z(J4) = ZERO
    DNM1 = Z(J4P2 + 2)
    DMIN = DNM1
    EMIN = ZERO
  else if (SAFMIN * Z(J4P2 + 2) < Z(J4 - 2) .and. &
 &         SAFMIN * Z(J4 - 2) < Z(J4P2 + 2)) then
    TEMP = Z(J4P2 + 2) / Z(J4 - 2)
    Z(J4) = Z(J4P2) * TEMP
    DNM1 = DNM2 * TEMP
  else
    Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
    DNM1 = Z(J4P2 + 2) * (DNM2 / Z(J4 - 2))
  end if
  DMIN = MIN(DMIN, DNM1)
!
  DMIN1 = DMIN
  J4 = J4 + 4
  J4P2 = J4 + 2 * PP - 1
  Z(J4 - 2) = DNM1 + Z(J4P2)
  if (Z(J4 - 2) == ZERO) then
    Z(J4) = ZERO
    DN = Z(J4P2 + 2)
    DMIN = DN
    EMIN = ZERO
  else if (SAFMIN * Z(J4P2 + 2) < Z(J4 - 2) .and. &
 &         SAFMIN * Z(J4 - 2) < Z(J4P2 + 2)) then
    TEMP = Z(J4P2 + 2) / Z(J4 - 2)
    Z(J4) = Z(J4P2) * TEMP
    DN = DNM1 * TEMP
  else
    Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
    DN = Z(J4P2 + 2) * (DNM1 / Z(J4 - 2))
  end if
  DMIN = MIN(DMIN, DN)
!
  Z(J4 + 2) = DN
  Z(4 * N0 - PP) = EMIN
  return
!
!     End of DLASQ6
!
end subroutine DLASQ6

