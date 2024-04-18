program main
  use mod_iolib
  use mod_progress_bar
  use mod_repeat_iterator
  use mod_cycle_iterator
  use mod_progbar_iterator
  use mod_forbar
  use mod_unittest
  implicit none
! type(unittest) :: u
  block
  type(repeat_iterator)  :: rep1, rep2
  type(cycle_iterator)   :: cyc
  type(progbar_iterator) :: prog
  type(repeat_iterator)  :: rep(3)
  type(forbar)           :: fbar1, fbar2, mbar
  integer                :: i
!
! call u%finish_and_terminate()
!
  rep1 = repeat_iterator(delay = 1000)
  rep2 = repeat_iterator(set=2, delay = 1000)
  cyc  = cycle_iterator(delay=1000)
  prog = progbar_iterator(limit=100000)
  fbar1 = forbar(cyc)
  fbar2 = forbar(prog)
  rep(1) = repeat_iterator(set=2, delay = 2500)
  rep(2) = repeat_iterator(set=1, delay = 1000)
  rep(3) = repeat_iterator(words=['*    ', '**   ', '***  ', ' *** ', '  ** ', '   * '], delay = 100)
  mbar = forbar(rep)
! do i = 1, 100000
!   call fbar2%update()
!   call fbar1%update()
!   write (*, '(A)', advance='NO') decorator(carret=.true.)
! end do
! print*
! print*
!
  do i = 1, 400000
    call mbar%update()
    write (*, '(A)', advance='NO') decorator(carret=.true.)
  end do
  print*
  print*
!
  end block
!contains
! subroutine test_forbar()
! end subroutine test_forbar
!
end program main
