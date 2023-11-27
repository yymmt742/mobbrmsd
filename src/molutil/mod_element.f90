module mod_element
  use mod_params, only: IK, RK, NL => NEWLINE
  implicit none
  private
  public :: element
  public :: element_formatted_string
  public :: ELEMENT_FORMATTED_CHARLEN
  public :: ELEMENT_FORMATTED_HEADER
!
  integer, parameter      :: element_formatted_charlen = 80
  character(*), parameter :: element_formatted_header  = &
 &   '  name  type  spec  atmnum            mass          charge'//NL
!
  type element
    sequence
    character(4) :: atmnam = 'X'
    character(4) :: atmtyp = 'X'
    character(4) :: nucspc = 'X'
    integer(IK)  :: atmnum
    real(RK)     :: atmass
    real(RK)     :: charge
  end type element
!
contains
!
  pure function element_formatted_string(this) result(res)
    type(element), intent(in)            :: this
    character(element_formatted_charlen) :: res
    write(res, '(2X,A,2X,A,2X,A,I8,2F16.9)') this
    res(element_formatted_charlen:element_formatted_charlen) = NL
  end function element_formatted_string
!
end module mod_element
