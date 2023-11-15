module mod_element
  use mod_params, only: IK, RK
  implicit none
  private
  public :: element
!
  type element
    sequence
    character(4) :: atmnam = 'X'
    character(4) :: atmtyp = 'X'
    character(4) :: nucspc = 'X'
    integer(IK)  :: atmnum = 0
    real(RK)     :: atmass = 1.0_RK
    real(RK)     :: charge = 0.0_RK
  end type element
!
end module mod_element
