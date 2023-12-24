program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_unittest
  use mod_tree
  implicit none
  type(unittest) :: u
!
  call u%init('test tree')
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    type(tree) :: t
!
    t = tree(1, 10, 4, [1,5,3,2])
!
  end subroutine test1
!
end program main
