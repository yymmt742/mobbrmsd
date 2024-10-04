!> \brief \b mobbrmsd_DLAS2 computes singular values of a 2-by-2 triangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DLAS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlas2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlas2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlas2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DLAS2( F, G, H, SSMIN, SSMAX )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   F, G, H, SSMAX, SSMIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DLAS2  computes the singular values of the 2-by-2 matrix
!>    [  F   G  ]
!>    [  0   H  ].
!> On return, SSMIN is the smaller singular value and SSMAX is the
!> larger singular value.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is DOUBLE PRECISION
!>          The (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is DOUBLE PRECISION
!>          The (1,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is DOUBLE PRECISION
!>          The (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] SSMIN
!> \verbatim
!>          SSMIN is DOUBLE PRECISION
!>          The smaller singular value.
!> \endverbatim
!>
!> \param[out] SSMAX
!> \verbatim
!>          SSMAX is DOUBLE PRECISION
!>          The larger singular value.
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
!> \ingroup OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Barring over/underflow, all output quantities are correct to within
!>  a few units in the last place (ulps), even in the absence of a guard
!>  digit in addition/subtraction.
!>
!>  In IEEE arithmetic, the code works correctly if one matrix element is
!>  infinite.
!>
!>  Overflow will not occur unless the largest singular value itself
!>  overflows, or is within a few ulps of overflow. (On machines with
!>  partial overflow, like the Cray, overflow may occur if the largest
!>  singular value is within a factor of 2 of overflow.)
!>
!>  Underflow is harmless if underflow is gradual. Otherwise, results
!>  may correspond to a matrix modified by perturbations of size near
!>  the underflow threshold.
!> \endverbatim
!>
!  =====================================================================
pure elemental subroutine mobbrmsd_DLAS2(F, G, H, SSMIN, SSMAX)
! use LA_CONSTANTS, only: RK => dp
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  real(RK), intent(in)  :: F, G, H
  real(RK), intent(out) :: SSMAX, SSMIN
!     ..
!
!  ====================================================================
!     ..
!     .. Local Scalars ..
  real(RK)             ::  AS, AT, AU, C, FA, FHMN, FHMX, GA, HA
!     ..
!     .. Intrinsic Functions ..
  intrinsic            :: ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     .. Parameters ..
! real(RK), parameter   :: ZERO = 0.0_RK
! real(RK), parameter   :: ONE = 1.0_RK
! real(RK), parameter   :: TWO = 2.0_RK
!
  FA = ABS(F)
  GA = ABS(G)
  HA = ABS(H)
  FHMN = MIN(FA, HA)
  FHMX = MAX(FA, HA)
  if (FHMN == ZERO) then
    SSMIN = ZERO
    if (FHMX == ZERO) then
      SSMAX = GA
    else
      SSMAX = MAX(FHMX, GA) * SQRT(ONE + (MIN(FHMX, GA) / MAX(FHMX, GA))**2)
    end if
  else
    if (GA < FHMX) then
      AS = ONE + FHMN / FHMX
      AT = (FHMX - FHMN) / FHMX
      AU = (GA / FHMX)**2
      C = TWO / (SQRT(AS * AS + AU) + SQRT(AT * AT + AU))
      SSMIN = FHMN * C
      SSMAX = FHMX / C
    else
      AU = FHMX / GA
      if (AU == ZERO) then
!
!              Avoid possible harmful underflow if exponent range
!              asymmetric (true SSMIN may not underflow even if
!              AU underflows)
!
        SSMIN = (FHMN * FHMX) / GA
        SSMAX = GA
      else
        AS = ONE + FHMN / FHMX
        AT = (FHMX - FHMN) / FHMX
        C = ONE / (SQRT(ONE + (AS * AU)**2) + SQRT(ONE + (AT * AU)**2))
        SSMIN = (FHMN * C) * AU
        SSMIN = SSMIN + SSMIN
        SSMAX = GA / (C + C)
      end if
    end if
  end if
  return
!
!     End of mobbrmsd_DLAS2
!
end subroutine mobbrmsd_DLAS2

