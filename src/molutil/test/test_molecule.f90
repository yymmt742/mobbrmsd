program main
  use mod_params, only: RK, IK, RZERO
  use mod_element
  use mod_element_set
  use mod_molecule
  use mod_molecular_topology
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call test1()
  call test2()
  call test3()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    type(element) :: e1, e2, e3
    print*,e1
    print*,e2
    print*,e3
  end subroutine test1
!
  subroutine test2()
    type(element_set) :: s1
    s1 = element_set()
    print*,s1
  end subroutine test2
!
  subroutine test3()
    type(molecule) :: s1
    s1 = molecule(3)
    print*,s1
    s1 = molecule(3, d=2, sym=[1, 3, 2])
    print*,s1
  end subroutine test3
!
end program main
