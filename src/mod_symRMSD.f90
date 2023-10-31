module mod_symRMSD
  use mod_unittest
  implicit none
!
contains
!
  subroutine X()
  type(unittest)               :: u
  logical,parameter            :: T=.TRUE., F=.FALSE.
  call u%init( 'test_unittest' )
  call u%assert(                T,           'assert               bool    0  ' )
  call u%finish_and_terminate()
  end subroutine
!
end module mod_symRMSD
