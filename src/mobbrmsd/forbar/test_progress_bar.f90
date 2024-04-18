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
  type(cycle_iterator)   :: cyc1, cyc2
  type(progbar_iterator) :: prog
  type(forbar)           :: fbar1, fbar2
  type(multi_forbar)     :: mbar
  integer                :: i
!
! call u%finish_and_terminate()
!
  rep1 = repeat_iterator(delay = 1000)
  rep2 = repeat_iterator(set=2, delay = 1000)
  cyc1 = cycle_iterator(word='***  ', delay = 100)
  cyc2 = cycle_iterator(delay=1000)
  prog = progbar_iterator(limit=100000, title='Test iterator   ', indent=4)
  fbar1 = forbar(cyc1)
  fbar2 = forbar(prog)
  mbar = multi_forbar()
  call mbar%add(prog)
  call mbar%add(rep2)

  i = 0
  do while(fbar1%running())
    write (*, '(A)', advance='NO') decorator(carret=.true.)
    i = i + 1
    if (i > 100000) exit
  end do
  print*
!
  i = 0
  do while(mbar%running())
    write (*, '(A)', advance='NO') decorator(carret=.true.)
    i = i + 1
  end do
  print*, i
!
  end block
!contains
! subroutine test_forbar()
! end subroutine test_forbar
!
end program main
