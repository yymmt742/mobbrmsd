!| mobbrmsd_LSAME returns .TRUE. if CA is the same letter as CB regardless of case.
!
!  Reference LSAME is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure elemental function mobbrmsd_LSAME(CA, CB)
  implicit none
!
!     .. Scalar Arguments ..
  character, intent(in) :: CA
!! CA specifies the single characters to be compared.
!!
  character, intent(in) :: CB
!! CB specifies the single characters to be compared.
!!
  logical               :: mobbrmsd_LSAME
!! Returns .TRUE. if CA is the same letter as CB regardless of case.
!!
  intrinsic :: ICHAR
  integer :: INTA, INTB, ZCODE
!
! Test if the characters are equal
!
  mobbrmsd_LSAME = CA == CB
  if (mobbrmsd_LSAME) return
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
      & INTA >= 145 .and. INTA <= 153 .or. &
      & INTA >= 162 .and. INTA <= 169) INTA = INTA + 64
    if (INTB >= 129 .and. INTB <= 137 .or. &
      & INTB >= 145 .and. INTB <= 153 .or. &
      & INTB >= 162 .and. INTB <= 169) INTB = INTB + 64
!
  else if (ZCODE == 218 .or. ZCODE == 250) then
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
    if (INTA >= 225 .and. INTA <= 250) INTA = INTA - 32
    if (INTB >= 225 .and. INTB <= 250) INTB = INTB - 32
  end if
  mobbrmsd_LSAME = INTA == INTB
!
  return
!
!     End of mobbrmsd_LSAME
!
end function mobbrmsd_LSAME

