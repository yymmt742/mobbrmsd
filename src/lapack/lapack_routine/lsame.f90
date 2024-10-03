!> \brief \b LSAME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION LSAME(CA,CB)
!
!       .. Scalar Arguments ..
!       CHARACTER CA,CB
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAME returns .TRUE. if CA is the same letter as CB regardless of
!> case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CA
!> \verbatim
!>          CA is CHARACTER!1
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CB is CHARACTER!1
!>          CA and CB specify the single characters to be compared.
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
!> \ingroup aux_blas
!
!  =====================================================================
pure elemental function LSAME(CA, CB)
  implicit none
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character, intent(in) :: CA, CB
  logical               :: LSAME
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
  intrinsic :: ICHAR
!     ..
!     .. Local Scalars ..
  integer :: INTA, INTB, ZCODE
!     ..
!
!     Test if the characters are equal
!
  LSAME = CA == CB
  if (LSAME) return
!
!     Now test for equivalence if both characters are alphabetic.
!
  ZCODE = ICHAR('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
  INTA = ICHAR(CA)
  INTB = ICHAR(CB)
!
  if (ZCODE == 90 .or. ZCODE == 122) then
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
    if (INTA >= 97 .and. INTA <= 122) INTA = INTA - 32
    if (INTB >= 97 .and. INTB <= 122) INTB = INTB - 32
!
  else if (ZCODE == 233 .or. ZCODE == 169) then
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
    if (INTA >= 129 .and. INTA <= 137 .or. &
&        INTA >= 145 .and. INTA <= 153 .or. &
&        INTA >= 162 .and. INTA <= 169) INTA = INTA + 64
    if (INTB >= 129 .and. INTB <= 137 .or. &
&        INTB >= 145 .and. INTB <= 153 .or. &
&        INTB >= 162 .and. INTB <= 169) INTB = INTB + 64
!
  else if (ZCODE == 218 .or. ZCODE == 250) then
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
    if (INTA >= 225 .and. INTA <= 250) INTA = INTA - 32
    if (INTB >= 225 .and. INTB <= 250) INTB = INTB - 32
  end if
  LSAME = INTA == INTB
!
  return
!
!     End of LSAME
!
end function LSAME
