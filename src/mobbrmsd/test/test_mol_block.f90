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
    print*,b%q
    print*,mol_block_nmol(b%q)
    print*,mol_block_napm(b%q)
    print*,mol_block_nsym(b%q)
!
    b = mol_block(8, 4, sym=RESHAPE([2, 4, 6, 8, 1, 3, 5, 7], [8, 1]))
    print'(*(I4))',b%q
    print*,mol_block_nmol(b%q)
    print*,mol_block_napm(b%q)
    print*,mol_block_nsym(b%q)
!
  end subroutine test1
!
end program main

