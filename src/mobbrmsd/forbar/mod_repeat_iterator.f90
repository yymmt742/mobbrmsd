module mod_repeat_iterator
  use mod_word_iterator
  implicit none
  private
  public :: repeat_iterator
!
!| NUMBER of ANIMATION frames
  integer, parameter  :: N_ANIMATION_A = 6
!
!| Character length of ANIMATION frame
  integer, parameter  :: L_ANIMATION_A = MAX( &
                       &                 LEN('__-¯¯'), &
                       &                 LEN('-__-¯'), &
                       &                 LEN('¯-__-'), &
                       &                 LEN('¯¯-__'), &
                       &                 LEN('-¯¯-_'), &
                       &                 LEN('_-¯¯-') &
                       &               )
!
!| ANIMATION strings
  character(L_ANIMATION_A), parameter :: ANIMATION_A1  = '__-¯¯'
  character(L_ANIMATION_A), parameter :: ANIMATION_A2  = '-__-¯'
  character(L_ANIMATION_A), parameter :: ANIMATION_A3  = '¯-__-'
  character(L_ANIMATION_A), parameter :: ANIMATION_A4  = '¯¯-__'
  character(L_ANIMATION_A), parameter :: ANIMATION_A5  = '-¯¯-_'
  character(L_ANIMATION_A), parameter :: ANIMATION_A6  = '_-¯¯-'
!
  character(L_ANIMATION_A), parameter :: ANIMATION_A(N_ANIMATION_A) = [ &
                                       &   ANIMATION_A1, &
                                       &   ANIMATION_A2, &
                                       &   ANIMATION_A3, &
                                       &   ANIMATION_A4, &
                                       &   ANIMATION_A5, &
                                       &   ANIMATION_A6  &
                                       &  ]
!
!| NUMBER of ANIMATION frames
  integer, parameter  :: N_ANIMATION_B = 10
!
!| Character length of ANIMATION frame
  integer, parameter  :: L_ANIMATION_B = MAX( &
                       &                 LEN('░░▒▒▓▓'), &
                       &                 LEN(' ░░▒▒▓'), &
                       &                 LEN('░ ░░▒▒'), &
                       &                 LEN('▒░ ░░▒'), &
                       &                 LEN('▒▒░ ░░'), &
                       &                 LEN('▓▒▒░ ░'), &
                       &                 LEN('▓▓▒▒░ '), &
                       &                 LEN('▒▓▓▒▒░'), &
                       &                 LEN('▒▒▓▓▒▒'), &
                       &                 LEN(' ▒▒▓▓▒') &
                       &               )
!
!| ANIMATION strings
  character(L_ANIMATION_B), parameter :: ANIMATION_B1  ='░░▒▒▓▓'
  character(L_ANIMATION_B), parameter :: ANIMATION_B2  =' ░░▒▒▓'
  character(L_ANIMATION_B), parameter :: ANIMATION_B3  ='░ ░░▒▒'
  character(L_ANIMATION_B), parameter :: ANIMATION_B4  ='▒░ ░░▒'
  character(L_ANIMATION_B), parameter :: ANIMATION_B5  ='▒▒░ ░░'
  character(L_ANIMATION_B), parameter :: ANIMATION_B6  ='▓▒▒░ ░'
  character(L_ANIMATION_B), parameter :: ANIMATION_B7  ='▓▓▒▒░ '
  character(L_ANIMATION_B), parameter :: ANIMATION_B8  ='▒▓▓▒▒░'
  character(L_ANIMATION_B), parameter :: ANIMATION_B9  ='▒▒▓▓▒▒'
  character(L_ANIMATION_B), parameter :: ANIMATION_B0  =' ▒▒▓▓▒'
!
  character(L_ANIMATION_B), parameter :: ANIMATION_B(N_ANIMATION_B) = [ &
                                     &   ANIMATION_B1, &
                                     &   ANIMATION_B2, &
                                     &   ANIMATION_B3, &
                                     &   ANIMATION_B4, &
                                     &   ANIMATION_B5, &
                                     &   ANIMATION_B6, &
                                     &   ANIMATION_B7, &
                                     &   ANIMATION_B8, &
                                     &   ANIMATION_B9, &
                                     &   ANIMATION_B0 &
                                     &  ]
!
!| Unicode BLOCK characters
!
  type, extends(word_iterator) :: repeat_iterator
    private
    integer                   :: counter
    integer                   :: wlen, tlen
    character(:), allocatable :: words
  contains
    procedure :: next => repeat_iterator_next
    final :: repeat_iterator_destroy
  end type repeat_iterator
!
  interface repeat_iterator
    module procedure repeat_iterator_new
  end interface repeat_iterator
!
contains
!
!| constructor
  pure function repeat_iterator_new(set, words) result(res)
    integer, intent(in), optional      :: set
    character(*), intent(in), optional :: words(:)
    type(repeat_iterator)              :: res
!
    res%counter = 0
!
    if (PRESENT(set)) then
      select case (set)
      case (1)
        res%wlen = L_ANIMATION_A
        res%tlen = N_ANIMATION_A * L_ANIMATION_A
        call alloc(res%wlen, res%tlen, ANIMATION_A, res%words, res%var)
        return
      case (2)
        res%wlen = L_ANIMATION_B
        res%tlen = N_ANIMATION_B * L_ANIMATION_B
        call alloc(res%wlen, res%tlen, ANIMATION_B, res%words, res%var)
        return
      end select
    else
      if (PRESENT(words)) then
        res%wlen = LEN(words)
        res%tlen = res%wlen * SIZE(words)
        call alloc(res%wlen, res%tlen, words, res%words, res%var)
        return
      end if
    end if
!
    res%wlen = L_ANIMATION_A
    res%tlen = N_ANIMATION_A * L_ANIMATION_A
    call alloc(res%wlen, res%tlen, ANIMATION_A, res%words, res%var)
!
  contains
    pure subroutine alloc(wlen, tlen, words, reswords, var)
      integer, intent(in)      :: wlen, tlen
      character(*), intent(in) :: words(:)
      character(:), allocatable, intent(inout) :: reswords
      character(:), allocatable, intent(inout) :: var
      integer                  :: i, j
        allocate(character(tlen) :: reswords)
        allocate(character(wlen) :: var)
        var = words(1)
        j = 0
        do i = 1, tlen, wlen
          j = j + 1
          reswords(i:i + wlen - 1) = words(j)
        end do
    end subroutine alloc
  end function repeat_iterator_new
!
!| iterator
  pure subroutine repeat_iterator_next(this)
    class(repeat_iterator), intent(inout) :: this
    integer                             :: i
      this%counter = MODULO(this%counter + this%wlen, this%tlen)
      do concurrent(i=1:this%wlen)
        block
          integer :: j
          j = this%counter + i
          this%var(i:i) = this%words(j:j)
        end block
      end do
  end subroutine repeat_iterator_next
!
!| destractor
  pure elemental subroutine repeat_iterator_destroy(this)
    type(repeat_iterator), intent(inout) :: this
    if (ALLOCATED(this%var)) deallocate (this%var)
    if (ALLOCATED(this%words)) deallocate (this%words)
  end subroutine repeat_iterator_destroy
!
end module mod_repeat_iterator

