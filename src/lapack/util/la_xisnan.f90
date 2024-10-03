# 1 "new/la_xisnan.F90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 399 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "new/la_xisnan.F90" 2
module LA_XISNAN
  implicit none
!
  private
  public :: LA_ISNAN, SISNAN, DISNAN
!
  interface LA_ISNAN
    module procedure SISNAN
    module procedure DISNAN
  end interface
!
contains
!
  pure elemental function SISNAN(x)
    use LA_CONSTANTS, only: SP
    real(SP), intent(in) :: x
    logical              :: SISNAN
    SISNAN = (x /= x)
  end function SISNAN
!
  pure elemental function DISNAN(x)
    use LA_CONSTANTS, only: DP
    real(DP), intent(in) :: x
    logical              :: DISNAN
    DISNAN = (x /= x)
  end function DISNAN
!
end module LA_XISNAN
