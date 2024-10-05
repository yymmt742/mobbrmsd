!| mobbrmsd_IEEECK is called from the mobbrmsd_ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Reference IEEECK is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- LAPACK auxiliary routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!     \date April 2012
!
pure elemental function mobbrmsd_IEEECK(ISPEC, ZERO, ONE)
  implicit none
  integer, intent(in)  :: ISPEC
!! SPecifies whether to test just for inifinity arithmetic
!! or whether to test for infinity and NaN arithmetic.
!!
!! = 0: Verify infinity arithmetic only.
!!
!! = 1: Verify infinity and NaN arithmetic.
!!
  real(RK), intent(in) ::  ZERO
!! Must contain the value 0.0.
!! This is passed to prevent the compiler from optimizing
!! away this code.
  real(RK), intent(in) :: ONE
!! Must contain the value 1.0.
!! This is passed to prevent the compiler from optimizing
!! away this code.
  integer              :: mobbrmsd_IEEECK
!! = 0:  Arithmetic failed to produce the correct answers
!!
!! = 1:  Arithmetic produced the correct answers
!!
  real(RK) :: NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, NEGZRO, NEWZRO, POSINF
!
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

