module mod_forbar
  use mod_optarg
  use mod_word_iterator
  use mod_iolib, only: &
 &      STDOUT => mod_iolib_STDOUT, &
 &      NULL_CHAR => mod_iolib_NULLCHR, &
 &      isatty, &
 &      decorate, &
 &      decorator
  implicit none
  private
  public :: forbar
!
  type forbar
    private
    class(word_iterator), allocatable :: iter
    class(word_iterator), allocatable :: multi(:)
  contains
    procedure :: update => forbar_update
    final :: forbar_destroy
  end type forbar
!
  interface forbar
    module procedure forbar_new, forbar_new_multi
  end interface forbar
!
contains
!
!| constructor
  pure function forbar_new(iter) result(res)
    class(word_iterator), intent(in) :: iter
    type(forbar)                     :: res
    ALLOCATE(res%iter, source=iter)
  end function forbar_new
!
!| constructor
  pure function forbar_new_multi(iter) result(res)
    class(word_iterator), intent(in) :: iter(:)
    type(forbar)                     :: res
    ALLOCATE(res%multi, source=iter)
  end function forbar_new_multi
!
  subroutine forbar_update(this)
    class(forbar), intent(inout) :: this
      if (ALLOCATED(this%iter)) then
        call this%iter%next()
        write (*, '(A)', advance='NO') this%iter%var
      elseif (ALLOCATED(this%multi)) then
        block
          integer :: i
          do concurrent(i=1:SIZE(this%multi))
            call this%multi(i)%next()
          end do
          do i = 1, SIZE(this%multi)
            write (*, '(A)', advance='NO') this%multi(i)%var
          end do
        end block
      end if
  end subroutine forbar_update
!
!| destractor
  pure elemental subroutine forbar_destroy(this)
    type(forbar), intent(inout) :: this
    if (ALLOCATED(this%iter)) deallocate (this%iter)
    if (ALLOCATED(this%multi)) deallocate (this%multi)
  end subroutine forbar_destroy
!
end module mod_forbar

