!> \brief \b SLARTG generates a plane rotation with real cosine and real sine.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARTG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARTG( F, G, CS, SN, R )
!
!       .. Scalar Arguments ..
!       REAL               CS, F, G, R, SN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARTG generate a plane rotation so that
!>
!>    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!>    [ -SN  CS  ]     [ G ]     [ 0 ]
!>
!> This is a slower, more accurate version of the BLAS1 routine SROTG,
!> with the following other differences:
!>    F and G are unchanged on return.
!>    If G=0, then CS=1 and SN=0.
!>    If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!>       floating point operations (saves work in SBDSQR when
!>       there are zeros on the diagonal).
!>
!> If F exceeds G in magnitude, CS will be positive.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is REAL
!>          The first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is REAL
!>          The second component of vector to be rotated.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is REAL
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is REAL
!>          The sine of the rotation.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is REAL
!>          The nonzero component of the rotated vector.
!>
!>  This version has a few statements commented out for thread safety
!>  (machine parameters are computed on each entry). 10 feb 03, SJH.
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
!> \date December 2016
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
pure elemental subroutine SLARTG(F, G, CS, SN, R)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
  real, intent(in)  :: F, G
  real, intent(out) :: CS, R, SN
!..
!
!  =====================================================================
!
!..Parameters..
  real, parameter :: ZERO = 0.0E0
  real, parameter :: ONE = 1.0E0
  real, parameter :: TWO = 2.0E0
!..
!..Local Scalars..
!logical FIRST
  integer :: CNT, I
  real :: EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCL
!..
  interface
! .. External Functions ..
    include 'slamch.h'
  end interface
!..
!..intrinsic Functions..
  intrinsic :: ABS, INT, LOG, MAX, SQRT
!..
!..save statement..
!save FIRST, SAFMX2, SAFMIN, SAFMN2
!..
!..data statements..
!data FIRST/.true./
!..
!..Executable Statements..
!
! if(FIRST) then
  SAFMIN = SLAMCH('S')
  EPS = SLAMCH('E')
  SAFMN2 = SLAMCH('B')**INT(LOG(SAFMIN / EPS) / LOG(SLAMCH('B')) / TWO)
  SAFMX2 = ONE / SAFMN2
! FIRST = .false.
! end if
  if (G == ZERO) then
    CS = ONE
    SN = ZERO
    R = F
  else if (F == ZERO) then
    CS = ZERO
    SN = ONE
    R = G
  else
    F1 = F
    G1 = G
    SCL = MAX(ABS(F1), ABS(G1))
    if (SCL >= SAFMX2) then
      CNT = 0
      do
        CNT = CNT + 1
        F1 = F1 * SAFMN2
        G1 = G1 * SAFMN2
        SCL = MAX(ABS(F1), ABS(G1))
        if (SCL < SAFMX2) exit
      end do
      R = SQRT(F1**2 + G1**2)
      CS = F1 / R
      SN = G1 / R
      do I = 1, CNT
        R = R * SAFMX2
      end do
    else if (SCL <= SAFMN2) then
      CNT = 0
      do
        CNT = CNT + 1
        F1 = F1 * SAFMX2
        G1 = G1 * SAFMX2
        SCL = MAX(ABS(F1), ABS(G1))
        if (SCL > SAFMN2) exit
      end do
      R = SQRT(F1**2 + G1**2)
      CS = F1 / R
      SN = G1 / R
      do I = 1, CNT
        R = R * SAFMN2
      end do
    else
      R = SQRT(F1**2 + G1**2)
      CS = F1 / R
      SN = G1 / R
    end if
    if (ABS(F) > ABS(G) .and. CS < ZERO) then
      CS = -CS
      SN = -SN
      R = -R
    end if
  end if
  return
  !
  !end of SLARTG
  !
end
