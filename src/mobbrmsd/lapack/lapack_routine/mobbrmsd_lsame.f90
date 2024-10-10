!| mobbrmsd_LSAME returns .TRUE. if CA is the same letter as CB regardless of case.
!
!  Reference LSAME is provided by [netlib](http://www.netlib.org/lapack/).
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  2024/10/06, Modified by yymmt742.
!
pure elemental function mobbrmsd_LSAME(CA, CB)
  implicit none
!
!     .. Scalar Arguments ..
  character, intent(in) :: CA
!! Single character to be compared.
!!
  character, intent(in) :: CB
!! Single character to be compared.
!!
  logical               :: mobbrmsd_LSAME
!! Returns TRUE if CA is the same letter as CB regardless of case.
!!
  character(*), parameter :: LET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
  integer, parameter      :: BW = ABS(IACHAR("a") - IACHAR("A"))
!
! Test if the characters are equal
!
  mobbrmsd_LSAME = CA == CB
  if (mobbrmsd_LSAME) return
  if (INDEX(LET, CA) <= 0) return
  mobbrmsd_LSAME = (ABS(ICHAR(CA) - ICHAR(CB)) == BW)
!
  return
!
! End of mobbrmsd_LSAME
!
end function mobbrmsd_LSAME

