!> \brief \b DLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARFG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgzfilename=/lapack/lapack_routine/dlarfg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zipfilename=/lapack/lapack_routine/dlarfg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txtfilename=/lapack/lapack_routine/dlarfg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       real(RK)           ::   ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       real(RK)           ::   X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARFG generates a real elementary reflector H of order n, such
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
!>          ALPHA is real(RK)           ::
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is real(RK)           :: array, dimension
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
!>          TAU is real(RK)           ::
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
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
pure subroutine DLARFG(N, ALPHA, X, INCX, TAU)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)     :: INCX, N
  real(RK), intent(inout) :: ALPHA
  real(RK), intent(out)   :: TAU
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: X(*)
!     ..
!
!  =====================================================================
!     .. Local Scalars ..
  integer              :: J, KNT
  real(RK)             :: BETA, RSAFMN, SAFMIN, XNORM
!
!     .. Intrinsic Functions ..
  intrinsic            :: ABS, SIGN
!     ..
!     .. Parameters ..
! real(RK), parameter   :: ZERO = 0.0_RK
! real(RK), parameter   :: ONE = 1.0_RK
!     ..
!     ..
! interface
!     .. External Functions ..
!   include 'dlamch.h'
!   include 'dlapy2.h'
!   include 'dnrm2.h'
!     .. External Subroutines ..
!   include 'dscal.h'
! end interface
!     ..
!     .. Executable Statements ..
!
  if (N <= 1) then
    TAU = ZERO
    return
  end if
!
  XNORM = DNRM2(N - 1, X, INCX)
!
  if (XNORM == ZERO) then
!
! H  =  I
!
    TAU = ZERO
  else
!
!        general case
!
    BETA = -SIGN(DLAPY2(ALPHA, XNORM), ALPHA)
    SAFMIN = DLAMCH('S') / DLAMCH('E')
    KNT = 0
    if (ABS(BETA) < SAFMIN) then
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
      RSAFMN = ONE / SAFMIN
10    continue
      KNT = KNT + 1
      call DSCAL(N - 1, RSAFMN, X, INCX)
      BETA = BETA * RSAFMN
      ALPHA = ALPHA * RSAFMN
      if ((ABS(BETA) < SAFMIN) .and. (KNT < 20)) GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
      XNORM = DNRM2(N - 1, X, INCX)
      BETA = -SIGN(DLAPY2(ALPHA, XNORM), ALPHA)
    end if
    TAU = (BETA - ALPHA) / BETA
    call DSCAL(N - 1, ONE / (ALPHA - BETA), X, INCX)
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
    do J = 1, KNT
      BETA = BETA * SAFMIN
    end do
    ALPHA = BETA
  end if
!
  return
!
!     End of DLARFG
!
end subroutine DLARFG
