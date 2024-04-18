module mod_timing_iterator
  use,intrinsic :: iso_fortran_env, only : INT64, REAL64
  use mod_optarg
  use mod_word_iterator
  implicit none
  private
  public :: timing_iterator
!
  character(*), parameter :: IFMT = '(A,I3,A,I2.2,A,I2.2,A,I3.3,A,F9.3,A)'
!
!| timing_iterator with ascii animation
  type, extends(word_iterator) :: timing_iterator
    private
    integer                   :: delay, idelay
    integer(INT64)            :: prev, start
    logical                   :: totty
    real(REAL64)              :: tosec
  contains
    procedure :: next    => timing_iterator_next
    procedure :: running => timing_iterator_running
    procedure :: to_tty  => timing_iterator_to_tty
    final :: timing_iterator_destroy
  end type timing_iterator
!
  interface timing_iterator
    module procedure timing_iterator_new
  end interface timing_iterator
!
contains
!
!| constructor
  pure function timing_iterator_new(delay, to_tty) result(res)
    integer, intent(in), optional :: delay
    logical, intent(in), optional :: to_tty
    type(timing_iterator)         :: res
    character(256)                :: var
!
    res%idelay = 0
    res%delay = optarg(delay, 1)
    res%totty = optarg(to_tty, .true.)

    res%prev  = -1
    write(var, IFMT) '[', 0,':', 0,':', 0,':', 0,', ', 0.0, ' iter/s]'
    allocate (res%var, source=TRIM(var))
!
  end function timing_iterator_new
!
!| iterator
  subroutine timing_iterator_next(this)
    class(timing_iterator), intent(inout) :: this
    integer(INT64)                       :: clock
    real(REAL64)                         :: time, iterpersec
    integer                              :: h, m, s, ms
!
      if (this%prev < 0) then
        call SYSTEM_CLOCK(this%start, clock)
        this%prev = this%start
        this%tosec = 1.0_REAL64 / clock
        return
      end if
!
      this%idelay = MODULO(this%idelay + 1, this%delay)
      if (this%idelay /= 0) return
!
      call SYSTEM_CLOCK(clock)
      if (clock > 1000 * 3600) clock = clock - 1000 * 3600
      time = real(clock - this%start, REAL64) * this%tosec
!
      s = INT(time)
      ms = 1000 * (time - s)
      h = s / 3600
      s = MODULO(s, 3600)
      m = s / 60
      s = MODULO(s, 60)
!
      iterpersec = this%delay / (real(clock - this%prev, REAL64) * this%tosec)
      this%prev = clock
!
      write (this%var, IFMT) '[', h, ':', m, ':', s, ':', ms, ', ', iterpersec, ' iter/s]'
!
  end subroutine timing_iterator_next
!
!| to_tty
  pure elemental function timing_iterator_to_tty(this) result(res)
    class(timing_iterator), intent(in) :: this
    logical                          :: res
    res = this%totty
  end function timing_iterator_to_tty
!
!| running
  pure elemental function timing_iterator_running(this) result(res)
    class(timing_iterator), intent(in) :: this
    logical                            :: res
    res = .true.
  end function timing_iterator_running
!
!| destractor
  pure elemental subroutine timing_iterator_destroy(this)
    type(timing_iterator), intent(inout) :: this
    if (ALLOCATED(this%var)) deallocate (this%var)
  end subroutine timing_iterator_destroy
!
end module mod_timing_iterator

