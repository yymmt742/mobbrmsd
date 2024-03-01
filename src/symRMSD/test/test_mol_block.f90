program main
  use mod_params, only: RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test mol_block')
  call test1()
  call u%finish_and_terminate()
!
contains
!
  subroutine test1()
    type(mol_block) :: b
!
    b = mol_block(10, 5)
    print*,mol_block_nmol(b)
    print*,mol_block_napm(b)
    print*,mol_block_nsym(b)
!
  end subroutine test1
!
end program main

