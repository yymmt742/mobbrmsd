program main
  use mod_params, only: D, setup_dimension, RK, IK, ONE => RONE, ZERO => RZERO
  use mod_mol_block
  use mod_lowerbound
  use mod_rotation_matrix
  use mod_testutil
  use mod_unittest
  implicit none
  type(unittest) :: u
!
  call u%init('test lowerbound')
!
  call setup_dimension(3)
  call test0()
!
  call u%finish_and_terminate()
!
contains
!
  subroutine test0()
    real(RK)        :: GC(10, 15), S(15), W(100), G, C(9)
    type(mol_block) :: b
    integer(IK)     :: i
!
    b = mol_block(8, 10)
    print*, lowerbound_worksize(b, 0)
    print*, lowerbound_worksize(b, 1)
    print*, lowerbound_worksize(b, 2)
    print*, lowerbound_worksize(b, 3)
!
    b = mol_block(8, 3)
    print*, lowerbound_worksize(b, 0)
    print*, lowerbound_worksize(b, 1)
    print*, lowerbound_worksize(b, 2)
    print*, lowerbound_worksize(b, 3)
!
    do i = 1, 15
      GC(:, i) = gcov(D, 8)
      call estimate_sdmin(GC(1, i), GC(2, i), W)
      S(i) = W(1)
    end do
!
    G = ZERO
    C = ZERO
!
    call lowerbound(0, b, G, C, S, W)
    print*,0, W(1)
    G = G + GC(1, 1)
    C = C + GC(2:, 1)
    call lowerbound(1, b, G, C, S([7,8,9,10,12,13,14,15]), W)
    print*,1, W(1)
    G = G + GC(1, 7)
    C = C + GC(2:, 7)
    call lowerbound(2, b, G, C, S([13,14,15]), W)
    print*,2, W(1)
    G = G + GC(1, 13)
    C = C + GC(2:, 13)
    call lowerbound(3, b, G, C, S, W)
    print*,3, W(1)
!
    G = ZERO
    C = ZERO
    G = G + GC(1, 3)
    C = C + GC(2:, 3)
    call lowerbound(1, b, G, C, S([6,7,9,10,11,12,14,15]), W)
    print*,1, W(1)
    G = G + GC(1, 9)
    C = C + GC(2:, 9)
    call lowerbound(2, b, G, C, S([11,12,15]), W)
    print*,2, W(1)
    G = G + GC(1, 15)
    C = C + GC(2:, 15)
    call lowerbound(3, b, G, C, S, W)
    print*,3, W(1)
!
  end subroutine test0
!
end program main
