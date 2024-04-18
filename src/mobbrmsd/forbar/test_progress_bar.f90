program main
  use mod_progress_bar
  use mod_repeat_iterator
  use mod_progbar_iterator
  use mod_unittest
  implicit none
  type(unittest) :: u
  type(repeat_iterator)  :: rep
  type(progbar_iterator) :: prog
! type(progress_bar) :: prog
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

  prog = progbar_iterator(limit=800)
  do i = 1, 810
    call prog%next()
    print'(A)', prog%var
  end do
  print*
!
end program main
