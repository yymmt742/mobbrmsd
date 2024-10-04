!> \brief \b mobbrmsd_SLASQ4 computes an approximation to the smallest eigenvalue using values of d from the previous transform. Used by sbdsqr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_SLASQ4 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq4.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq4.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq4.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_SLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN,
!                          DN1, DN2, TAU, TTYPE, G )
!
!       .. Scalar Arguments ..
!       INTEGER            I0, N0, N0IN, PP, TTYPE
!       REAL               DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
!       ..
!       .. Array Arguments ..
!       REAL               Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_SLASQ4 computes an approximation TAU to the smallest eigenvalue
!> using values of d from the previous transform.
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
!>          Z is REAL array, dimension ( 4*N0 )
!>        Z holds the qd array.
!> \endverbatim
!>
!> \param[in] PP
!> \verbatim
!>          PP is INTEGER
!>        PP=0 for ping, PP=1 for pong.
!> \endverbatim
!>
!> \param[in] N0IN
!> \verbatim
!>          N0IN is INTEGER
!>        The value of N0 at start of EIGTEST.
!> \endverbatim
!>
!> \param[in] DMIN
!> \verbatim
!>          DMIN is REAL
!>        Minimum value of d.
!> \endverbatim
!>
!> \param[in] DMIN1
!> \verbatim
!>          DMIN1 is REAL
!>        Minimum value of d, excluding D( N0 ).
!> \endverbatim
!>
!> \param[in] DMIN2
!> \verbatim
!>          DMIN2 is REAL
!>        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!> \endverbatim
!>
!> \param[in] DN
!> \verbatim
!>          DN is REAL
!>        d(N)
!> \endverbatim
!>
!> \param[in] DN1
!> \verbatim
!>          DN1 is REAL
!>        d(N-1)
!> \endverbatim
!>
!> \param[in] DN2
!> \verbatim
!>          DN2 is REAL
!>        d(N-2)
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL
!>        This is the shift.
!> \endverbatim
!>
!> \param[out] TTYPE
!> \verbatim
!>          TTYPE is INTEGER
!>        Shift type.
!> \endverbatim
!>
!> \param[in,out] G
!> \verbatim
!>          G is REAL
!>        G is passed as an argument in order to save its value between
!>        calls to mobbrmsd_SLASQ4.
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
!> \date June 2016
!
!> \ingroup auxOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  CNST1 = 9/16
!> \endverbatim
!>
!  =====================================================================
pure subroutine mobbrmsd_SLASQ4(I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, &
               &       DN1, DN2, TAU, TTYPE, G)
  implicit none
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
  integer, intent(in)  :: I0, N0, N0IN, PP
  integer, intent(out) :: TTYPE
  real(RK), intent(in)     :: DMIN, DMIN1, DMIN2, DN, DN1, DN2
  real(RK), intent(inout)  :: G
  real(RK), intent(out)    :: TAU
!..
!..Array Arguments..
  real(RK), intent(in)     :: Z(*)
!..
!
!  =====================================================================
!..
!..Local Scalars..
  integer :: I4, NN, NP
  real(RK) :: A2, B1, B2, GAM, GAP1, GAP2, S
!
!..Parameters..
  real(RK), parameter :: CNST1 = 0.5630_RK
  real(RK), parameter :: CNST2 = 1.010_RK
  real(RK), parameter :: CNST3 = 1.050_RK
! real(RK), parameter :: ZERO = 0.0E0
! real(RK), parameter :: QURTR = 0.250E0
! real(RK), parameter :: HALF = 0.5E0
! real(RK), parameter :: THIRD = 0.3330E0
! real(RK), parameter :: ONE = 1.0E+0
! real(RK), parameter :: TWO = 2.0E0
! real(RK), parameter :: HUNDRD = 100.0E0
!..
!..intrinsic Functions..
  intrinsic :: MAX, MIN, SQRT
!..
!..Executable Statements..
!
!A negative DMIN forces the shift to take that absolute value
!TTYPE records the type of shift.
!
  if (DMIN <= ZERO) then
    TAU = -DMIN
    TTYPE = -1
    return
  end if
!
  NN = 4 * N0 + PP
  if (N0IN == N0) then
    !
    !No eigenvalues deflated.
    !
    if (DMIN == DN .or. DMIN == DN1) then
      !
      B1 = SQRT(Z(NN - 3)) * SQRT(Z(NN - 5))
      B2 = SQRT(Z(NN - 7)) * SQRT(Z(NN - 9))
      A2 = Z(NN - 7) + Z(NN - 5)
      !
      !Cases 2 and 3.
      !
      if (DMIN == DN .and. DMIN1 == DN1) then
        GAP2 = DMIN2 - A2 - DMIN2 * QURTR
        if (GAP2 > ZERO .and. GAP2 > B2) then
          GAP1 = A2 - DN - (B2 / GAP2) * B2
        else
          GAP1 = A2 - DN - (B1 + B2)
        end if
        if (GAP1 > ZERO .and. GAP1 > B1) then
          S = MAX(DN - (B1 / GAP1) * B1, HALF * DMIN)
          TTYPE = -2
        else
          S = ZERO
          if (DN > B1) S = DN - B1
          if (A2 > (B1 + B2)) S = MIN(S, A2 - (B1 + B2))
          S = MAX(S, THIRD * DMIN)
          TTYPE = -3
        end if
      else
        !
        !case 4.
        !
        TTYPE = -4
        S = QURTR * DMIN
        if (DMIN == DN) then
          GAM = DN
          A2 = ZERO
          if (Z(NN - 5) > Z(NN - 7)) return
          B2 = Z(NN - 5) / Z(NN - 7)
          NP = NN - 9
        else
          NP = NN - 2 * PP
          GAM = DN1
          if (Z(NP - 4) > Z(NP - 2)) return
          A2 = Z(NP - 4) / Z(NP - 2)
          if (Z(NN - 9) > Z(NN - 11)) return
          B2 = Z(NN - 9) / Z(NN - 11)
          NP = NN - 13
        end if
        !
        !Approximate contribution to norm squared from I < NN - 1.
        !
        A2 = A2 + B2
        do I4 = NP, 4 * I0 - 1 + PP, -4
          if (B2 == ZERO) exit
          B1 = B2
          if (Z(I4) > Z(I4 - 2)) return
          B2 = B2 * (Z(I4) / Z(I4 - 2))
          A2 = A2 + B2
          if (HUNDRD * MAX(B2, B1) < A2 .or. CNST1 < A2) exit
        end do
        A2 = CNST3 * A2
        !
        !Rayleigh quotient residual bound.
        !
        if (A2 < CNST1) S = GAM * (ONE - SQRT(A2)) / (ONE + A2)
      end if
    else if (DMIN == DN2) then
      !
      !case 5.
      !
      TTYPE = -5
      S = QURTR * DMIN
      !
      !Compute contribution to norm squared from I > NN - 2.
      !
      NP = NN - 2 * PP
      B1 = Z(NP - 2)
      B2 = Z(NP - 6)
      GAM = DN2
      if (Z(NP - 8) > B2 .or. Z(NP - 4) > B1) return
      A2 = (Z(NP - 8) / B2) * (ONE + Z(NP - 4) / B1)
      !
      !Approximate contribution to norm squared from I < NN - 2.
      !
      if (N0 - I0 > 2) then
        B2 = Z(NN - 13) / Z(NN - 15)
        A2 = A2 + B2
        do I4 = NN - 17, 4 * I0 - 1 + PP, -4
          if (B2 == ZERO) exit
          B1 = B2
          if (Z(I4) > Z(I4 - 2)) return
          B2 = B2 * (Z(I4) / Z(I4 - 2))
          A2 = A2 + B2
          if (HUNDRD * MAX(B2, B1) < A2 .or. CNST1 < A2) exit
        end do
        A2 = CNST3 * A2
      end if
      !
      if (A2 < CNST1) S = GAM * (ONE - SQRT(A2)) / (ONE + A2)
    else
      !
      !case 6, no information to guide us.
      !
      if (TTYPE == -6) then
        G = G + THIRD * (ONE - G)
      else if (TTYPE == -18) then
        G = QURTR * THIRD
      else
        G = QURTR
      end if
      S = G * DMIN
      TTYPE = -6
    end if
    !
  else if (N0IN == (N0 + 1)) then
    !
    !One eigenvalue just deflated.use DMIN1, DN1 for DMIN and DN.
    !
    if (DMIN1 == DN1 .and. DMIN2 == DN2) then
      !
      !Cases 7 and 8.
      !
      TTYPE = -7
      S = THIRD * DMIN1
      if (Z(NN - 5) > Z(NN - 7)) return
      B1 = Z(NN - 5) / Z(NN - 7)
      B2 = B1
      if (B2 /= ZERO) then
        do I4 = 4 * N0 - 9 + PP, 4 * I0 - 1 + PP, -4
          A2 = B1
          if (Z(I4) > Z(I4 - 2)) return
          B1 = B1 * (Z(I4) / Z(I4 - 2))
          B2 = B2 + B1
          if (HUNDRD * MAX(B1, A2) < B2) exit
        end do
      end if
      B2 = SQRT(CNST3 * B2)
      A2 = DMIN1 / (ONE + B2**2)
      GAP2 = HALF * DMIN2 - A2
      if (GAP2 > ZERO .and. GAP2 > B2 * A2) then
        S = MAX(S, A2 * (ONE - CNST2 * A2 * (B2 / GAP2) * B2))
      else
        S = MAX(S, A2 * (ONE - CNST2 * B2))
        TTYPE = -8
      end if
    else
      !
      !case 9.
      !
      S = QURTR * DMIN1
      if (DMIN1 == DN1) S = HALF * DMIN1
      TTYPE = -9
    end if
    !
  else if (N0IN == (N0 + 2)) then
    !
    !Two eigenvalues deflated.use DMIN2, DN2 for DMIN and DN.
    !
    !Cases 10 and 11.
    !
    if (DMIN2 == DN2 .and. TWO * Z(NN - 5) < Z(NN - 7)) then
      TTYPE = -10
      S = THIRD * DMIN2
      if (Z(NN - 5) > Z(NN - 7)) return
      B1 = Z(NN - 5) / Z(NN - 7)
      B2 = B1
      if (B2 /= ZERO) then
        do I4 = 4 * N0 - 9 + PP, 4 * I0 - 1 + PP, -4
          if (Z(I4) > Z(I4 - 2)) return
          B1 = B1 * (Z(I4) / Z(I4 - 2))
          B2 = B2 + B1
          if (HUNDRD * B1 < B2) exit
        end do
      end if
      B2 = SQRT(CNST3 * B2)
      A2 = DMIN2 / (ONE + B2**2)
      GAP2 = Z(NN - 7) + Z(NN - 9) - SQRT(Z(NN - 11)) * SQRT(Z(NN - 9)) - A2
      if (GAP2 > ZERO .and. GAP2 > B2 * A2) then
        S = MAX(S, A2 * (ONE - CNST2 * A2 * (B2 / GAP2) * B2))
      else
        S = MAX(S, A2 * (ONE - CNST2 * B2))
      end if
    else
      S = QURTR * DMIN2
      TTYPE = -11
    end if
  else if (N0IN > (N0 + 2)) then
    !
    !case 12, more than two eigenvalues deflated.No information.
    !
    S = ZERO
    TTYPE = -12
  end if
  !
  TAU = S
  return
  !
  !end of mobbrmsd_SLASQ4
  !
end
