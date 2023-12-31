program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_optarg
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test optarg logical')
  call test1(.FALSE.)
  call test1(.TRUE., .TRUE.)
  call test1(.FALSE., .FALSE.)
!
  call test2(.TRUE., .TRUE.)
  call test2(.FALSE., .FALSE.)
  call test2(.TRUE., .TRUE., .TRUE.)
  call test2(.TRUE., .FALSE., .TRUE.)
  call test2(.FALSE., .TRUE., .FALSE.)
  call test2(.FALSE., .FALSE., .FALSE.)
!
  call u%init('test optarg integer')
!
  call test3(1, 1)
  call test3(2, 1, 2)
!
  call u%init('test optarg real')
!
  call test4(1._RK, 1._RK)
  call test4(2._RK, 1._RK, 2._RK)
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1(T, l)
    logical, intent(in)           :: T
    logical, intent(in), optional :: l
!
    call u%assert_equal(T, optarg(l), 'optarg_l0 ')
!
  end subroutine test1
!
  subroutine test2(T, def, l)
    logical, intent(in)           :: T, def
    logical, intent(in), optional :: l
!
    call u%assert_equal(T, optarg(l, def), 'optarg_l0_')
!
  end subroutine test2
!
  subroutine test3(T, def, l)
    integer(IK), intent(in)           :: T, def
    integer(IK), intent(in), optional :: l
!
    call u%assert_equal(T, optarg(l, def), 'optarg_i0')
!
  end subroutine test3
!
  subroutine test4(T, def, l)
    real(RK), intent(in)           :: T, def
    real(RK), intent(in), optional :: l
!
    call u%assert_almost_equal(T, optarg(l, def), 'optarg_r0')
!
  end subroutine test4
!
end program main
