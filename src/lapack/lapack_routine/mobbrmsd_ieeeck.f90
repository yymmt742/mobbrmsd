!> \brief \b mobbrmsd_IEEECK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
! http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_IEEECK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!   INTEGER          FUNCTION mobbrmsd_IEEECK( ISPEC, ZERO, ONE )
!
!   .. Scalar Arguments ..
!   INTEGER            ISPEC
!   REAL               ONE, ZERO
!   ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_IEEECK is called from the mobbrmsd_ILAENV to verify that Infinity and
!> possibly NaN arithmetic is safe (i.e. will not trap).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          SPecifies whether to test just for inifinity arithmetic
!>          or whether to test for infinity and NaN arithmetic.
!>          = 0: Verify infinity arithmetic only.
!>          = 1: Verify infinity and NaN arithmetic.
!> \endverbatim
!>
!> \param[in] ZERO
!> \verbatim
!>          ZERO is REAL
!>          Must contain the value 0.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!> \endverbatim
!>
!> \param[in] ONE
!> \verbatim
!>          ONE is REAL
!>          Must contain the value 1.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!>
!>  RETURN VALUE:  INTEGER
!>          = 0:  Arithmetic failed to produce the correct answers
!>          = 1:  Arithmetic produced the correct answers
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
!  =====================================================================
pure elemental function mobbrmsd_IEEECK(ISPEC, ZERO, ONE)
! use LA_CONSTANTS, only: SP
  implicit none
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
! .. Scalar Arguments ..
  integer, intent(in)  :: ISPEC
  real(RK), intent(in) :: ONE, ZERO
  integer              :: mobbrmsd_IEEECK
! ..
!
!  =====================================================================
!
! .. Local Scalars ..
  real(RK) :: NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, NEGZRO, NEWZRO, POSINF
! ..
! .. Executable Statements ..
  mobbrmsd_IEEECK = 1
!
  POSINF = ONE / ZERO
  if (POSINF <= ONE) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  NEGINF = -ONE / ZERO
  if (NEGINF >= ZERO) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  NEGZRO = ONE / (NEGINF + ONE)
  if (NEGZRO /= ZERO) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  NEGINF = ONE / NEGZRO
  if (NEGINF >= ZERO) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  NEWZRO = NEGZRO + ZERO
  if (NEWZRO /= ZERO) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  POSINF = ONE / NEWZRO
  if (POSINF <= ONE) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  NEGINF = NEGINF * POSINF
  if (NEGINF >= ZERO) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  POSINF = POSINF * POSINF
  if (POSINF <= ONE) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
! Return if we were only asked to check infinity arithmetic
!
  if (ISPEC == 0) return
!
  NAN1 = POSINF + NEGINF
!
  NAN2 = POSINF / NEGINF
!
  NAN3 = POSINF / POSINF
!
  NAN4 = POSINF * ZERO
!
  NAN5 = NEGINF * NEGZRO
!
  NAN6 = NAN5 * ZERO
!
  if (NAN1 == NAN1) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  if (NAN2 == NAN2) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  if (NAN3 == NAN3) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  if (NAN4 == NAN4) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  if (NAN5 == NAN5) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  if (NAN6 == NAN6) then
    mobbrmsd_IEEECK = 0
    return
  end if
!
  return
end function mobbrmsd_IEEECK
