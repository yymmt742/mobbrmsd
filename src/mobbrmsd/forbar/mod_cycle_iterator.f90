module mod_cycle_iterator
  use mod_optarg
  use mod_word_iterator
  implicit none
  private
  public :: cycle_iterator
!
  character(*), parameter :: ANIMATION_A = '_-¯¯-_'
  character(*), parameter :: ANIMATION_B = '------>  '
!
!| cycle_iterator with ascii animation
  type, extends(word_iterator) :: cycle_iterator
    private
    integer                   :: delay, idelay
  contains
    procedure :: next => cycle_iterator_next
    final :: cycle_iterator_destroy
  end type cycle_iterator
!
  interface cycle_iterator
    module procedure cycle_iterator_new
  end interface cycle_iterator
!
contains
!
!| constructor
  pure function cycle_iterator_new(set, word, delay) result(res)
    integer, intent(in), optional      :: set
    character(*), intent(in), optional :: word
    integer, intent(in), optional      :: delay
    type(cycle_iterator)               :: res
!
    res%idelay = 0
    res%delay = optarg(delay, 1)
!
    if (PRESENT(set)) then
      select case (set)
      case (1)
        allocate(res%var, source=ANIMATION_A)
      case (2)
        allocate(res%var, source=ANIMATION_B)
      end select
    else
      if (PRESENT(word)) then
        allocate(res%var, source=word)
      else
        allocate(res%var, source=ANIMATION_A)
      end if
    end if
  end function cycle_iterator_new
!
!| iterator
  pure subroutine cycle_iterator_next(this)
    class(cycle_iterator), intent(inout) :: this
    character                            :: tmp
    integer                              :: i, l
      this%idelay = MODULO(this%idelay + 1, this%delay)
      if (this%idelay /= 0) return
      l = LEN(this%var)
      tmp = this%var(l:l)
      do i = l, 2, -1
        this%var(i:i) = this%var(i - 1:i - 1)
      end do
      this%var(1:1) = tmp
  end subroutine cycle_iterator_next
!
!| destractor
  pure elemental subroutine cycle_iterator_destroy(this)
    type(cycle_iterator), intent(inout) :: this
    if (ALLOCATED(this%var)) deallocate (this%var)
  end subroutine cycle_iterator_destroy
!
end module mod_cycle_iterator

