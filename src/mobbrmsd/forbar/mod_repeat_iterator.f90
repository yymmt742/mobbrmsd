module mod_repeat_iterator
  use mod_optarg
  use mod_word_iterator
  implicit none
  private
  public :: repeat_iterator
!
!| NUMBER of ANIMATION frames
  integer, parameter  :: N_ANIMATION_A = 50
!
!| Character length of ANIMATION frame
  integer, parameter  :: L_ANIMATION_A = MAX( &
                       &                 LEN('▁▁'), &
                       &                 LEN('▂▂'), &
                       &                 LEN('▃▃'), &
                       &                 LEN('▄▄'), &
                       &                 LEN('▅▅'), &
                       &                 LEN('▆▆'), &
                       &                 LEN('▇▇'), &
                       &                 LEN('██'), &
                       &                 LEN('█▉'), &
                       &                 LEN('█▊'), &
                       &                 LEN('█▋'), &
                       &                 LEN('█▌'), &
                       &                 LEN('█▍'), &
                       &                 LEN('█▎'), &
                       &                 LEN('█▏') &
                       &               )
!
  character(L_ANIMATION_A), parameter :: A0 = '▁▁'
  character(L_ANIMATION_A), parameter :: A1 = '▂▂'
  character(L_ANIMATION_A), parameter :: A2 = '▃▃'
  character(L_ANIMATION_A), parameter :: A3 = '▄▄'
  character(L_ANIMATION_A), parameter :: A4 = '▅▅'
  character(L_ANIMATION_A), parameter :: A5 = '▆▆'
  character(L_ANIMATION_A), parameter :: A6 = '▇▇'
  character(L_ANIMATION_A), parameter :: A7 = '██'
  character(L_ANIMATION_A), parameter :: A8 = '█▉'
  character(L_ANIMATION_A), parameter :: A9 = '█▊'
  character(L_ANIMATION_A), parameter :: AA = '█▋'
  character(L_ANIMATION_A), parameter :: AB = '█▌'
  character(L_ANIMATION_A), parameter :: AC = '█▍'
  character(L_ANIMATION_A), parameter :: AD = '█▎'
  character(L_ANIMATION_A), parameter :: AE = '█▏'
!
!| ANIMATION strings
  character(L_ANIMATION_A), parameter :: ANIMATION_A(N_ANIMATION_A) = [ &
                                                      &  A0,  &
                                                      &  A0,  &
                                                      &  A0,  &
                                                      &  A0,  &
                                                      &  A0,  &
                                                      &  A0,  &
                                                      &  A0,  &
                                                      &  A0,  &
                                                      &  A0,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A2,  &
                                                      &  A3,  &
                                                      &  A4,  &
                                                      &  A5,  &
                                                      &  A6,  &
                                                      &  A7,  &
                                                      &  A8,  &
                                                      &  A9,  &
                                                      &  AA,  &
                                                      &  AB,  &
                                                      &  AC,  &
                                                      &  AD,  &
                                                      &  AE,  &
                                                      &  AD,  &
                                                      &  AC,  &
                                                      &  AB,  &
                                                      &  AA,  &
                                                      &  A9,  &
                                                      &  A8,  &
                                                      &  A7,  &
                                                      &  A6,  &
                                                      &  A5,  &
                                                      &  A4,  &
                                                      &  A3,  &
                                                      &  A2,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1,  &
                                                      &  A1   &
                                                      &  ]

!| NUMBER of ANIMATION frames
  integer, parameter  :: N_ANIMATION_B = 28
!
!| Character length of ANIMATION frame
  integer, parameter  :: L_ANIMATION_B = MAX( &
                       &                 LEN('...    '), &
                       &                 LEN('.. .   '), &
                       &                 LEN('. . .  '), &
                       &                 LEN('. .  . '), &
                       &                 LEN('.  .  .'), &
                       &                 LEN(' .  . .'), &
                       &                 LEN('   . ..'), &
                       &                 LEN('    ...') &
                       &               )
!
!| ANIMATION strings
!
  character(L_ANIMATION_B), parameter :: ANIMATION_B(N_ANIMATION_B) = [ &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '...    ', &
                                     &   '.. .   ', &
                                     &   '.. .   ', &
                                     &   '.. .   ', &
                                     &   '.. .   ', &
                                     &   '.. .   ', &
                                     &   '.. .   ', &
                                     &   '. . .  ', &
                                     &   '. .  . ', &
                                     &   '.   . .', &
                                     &   '  .  ..', &
                                     &   '   . ..', &
                                     &   '    ...', &
                                     &   '  . .. ', &
                                     &   '.. .   '  &
                                     &  ]
!
!| repeat_iterator with ascii animation
  type, extends(word_iterator) :: repeat_iterator
    private
    integer                   :: counter
    integer                   :: delay, idelay
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
  pure function repeat_iterator_new(set, words, delay) result(res)
    integer, intent(in), optional      :: set
    character(*), intent(in), optional :: words(:)
    integer, intent(in), optional      :: delay
    type(repeat_iterator)              :: res
!
    res%counter = 0
    res%idelay = 0
    res%delay = optarg(delay, 1)
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
!
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
    integer                               :: i
      this%idelay = MODULO(this%idelay + 1, this%delay)
      if (this%idelay /= 0) return
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

