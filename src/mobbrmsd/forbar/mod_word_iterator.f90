module mod_word_iterator
  implicit none
  private
  public :: word_iterator
!
  type word_iterator
    private
    character(:), allocatable, public :: var
  contains
    procedure :: next => word_iterator_next
    final :: word_iterator_destroy
  end type word_iterator
!
  interface word_iterator
    module procedure :: word_iterator_new
  end interface word_iterator
!
contains
!
!| constructor
  pure function word_iterator_new() result(res)
    type(word_iterator) :: res
    allocate(character(0) :: res%var)
  end function word_iterator_new
!
!| iterator
  pure subroutine word_iterator_next(this)
    class(word_iterator), intent(inout) :: this
  end subroutine word_iterator_next
!
!| destractor
  pure elemental subroutine word_iterator_destroy(this)
    type(word_iterator), intent(inout) :: this
    if (ALLOCATED(this%var)) deallocate (this%var)
  end subroutine word_iterator_destroy
!
end module mod_word_iterator

