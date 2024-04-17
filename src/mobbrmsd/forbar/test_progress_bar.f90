program main
  use mod_progress_bar
  use mod_repeat_iterator
  use mod_unittest
  implicit none
  type(unittest) :: u
  type(repeat_iterator) :: rep
  type(progress_bar) :: prog
  integer        :: i
!
  call u%finish_and_terminate()
!
  rep = repeat_iterator()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
!
  rep = repeat_iterator(set=2)
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var
  call rep%next()
  print*, rep%var

  prog = progress_bar(limit=110000, title='Title', indent=4)
  do i = 1, 1100000001
    call prog%update(i/10000)
  end do
  print*
!
end program main
