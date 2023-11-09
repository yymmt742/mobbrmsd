program main
  use mod_params, only: RK, IK, RZERO
  use mod_rmsd
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call test1()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
!
    real(RK) :: X(12)
    real(RK) :: Y(12)
!
    call u%init('test rmsd')
!
    call random_number(X)
    call random_number(Y)
    print *, rmsd(3, 4, X, Y)
!
  end subroutine test1
!
end program main
