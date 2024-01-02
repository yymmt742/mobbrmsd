module driver
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_symRMSD
  use mod_mol_block
  implicit none
  private
  public hello
!
  type(symRMSD_input) :: inp
  type(symRMSD)       :: sym
!
contains
!
  subroutine hello()
    print'(A)', 'Hello'
  end subroutine hello
!
! subroutine load(arg, ios)
!   integer(IK)              :: unit, nn(2)
!
!   open (NEWUNIT=unit, FILE=arg, STATUS='OLD', IOSTAT=ios)
!   if (ios > 0) stop
!   read (unit, *, IOSTAT=ios) d
!   if (ios > 0) stop
!   close (unit, IOSTAT=ios)
!
!   call ip%copy(d%ip)
!   n = d%nstate()
!   nn = [n, n]
!   ip = d%ip%index_pointer()
!
! end subroutine load
!
end module driver
