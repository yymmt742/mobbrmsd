!> \brief \b SLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARFG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL               ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       REAL               X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARFG generates a real elementary reflector H of order n, such
!> that
!>
!>       H * ( alpha ) = ( beta ),   H**T * H = I.
!>           (   x   )   (   0  )
!>
!> where alpha and beta are scalars, and x is an (n-1)-element real
!> vector. H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
!>                     ( v )
!>
!> where tau is a real scalar and v is a real (n-1)-element
!> vector.
!>
!> If the elements of x are all zero, then tau = 0 and H is taken to be
!> the unit matrix.
!>
!> Otherwise  1 <= tau <= 2.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the elementary reflector.
!> \endverbatim
!>
!> \param[in,out] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is REAL array, dimension
!>                         (1+(N-2)*abs(INCX))
!>          On entry, the vector x.
!>          On exit, it is overwritten with the vector v.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL
!>          The value tau.
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
!> \date November 2017
!
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
pure subroutine SLARFG(N, ALPHA, X, INCX, TAU)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
  integer, intent(in) :: INCX, N
  real, intent(inout) :: ALPHA
  real, intent(out)   :: TAU
!..
!..Array Arguments..
  real, intent(inout) :: X(*)
!..
!
!  =====================================================================
!
!..Parameters..
  real, parameter :: ZERO = 0.0E0
  real, parameter :: ONE = 1.0E0
!..
!..Local Scalars..
  integer :: J, KNT
  real :: BETA, RSAFMN, SAFMIN, XNORM
!..
!..intrinsic Functions..
  intrinsic :: ABS, SIGN
!..
  interface
!..external Functions..
    include 'slamch.h'
    include 'slapy2.h'
    include 'snrm2.h'
!..
!..external Subroutines..
    include 'sscal.h'
  end interface
!..
!..Executable Statements..
!
  if (N <= 1) then
    TAU = ZERO
    return
  end if
!
  XNORM = SNRM2(N - 1, X, INCX)
!
  if (XNORM == ZERO) then
!
! H = I
!
    TAU = ZERO
  else
!
! general case
!
    BETA = -SIGN(SLAPY2(ALPHA, XNORM), ALPHA)
    SAFMIN = SLAMCH('S') / SLAMCH('E')
    KNT = 0
    if (ABS(BETA) < SAFMIN) then
!
! XNORM, BETA may be inaccurate; scale X and recompute them
!
      RSAFMN = ONE / SAFMIN
      do
!10    continue
        KNT = KNT + 1
        call SSCAL(N - 1, RSAFMN, X, INCX)
        BETA = BETA * RSAFMN
        ALPHA = ALPHA * RSAFMN
!       if ((ABS(BETA) < SAFMIN) .and. (KNT < 20)) GO TO 10
        if ((ABS(BETA) >= SAFMIN) .or. (KNT >= 20)) exit
      end do
!
! New BETA is at most 1, at least SAFMIN
!
      XNORM = SNRM2(N - 1, X, INCX)
      BETA = -SIGN(SLAPY2(ALPHA, XNORM), ALPHA)
    end if
    TAU = (BETA - ALPHA) / BETA
    call SSCAL(N - 1, ONE / (ALPHA - BETA), X, INCX)
!
! if ALPHA is subnormal, it may lose relative accuracy
!
    do J = 1, KNT
      BETA = BETA * SAFMIN
    end do
    ALPHA = BETA
  end if
!
  return
!
! end of SLARFG
!
end
