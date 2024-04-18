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
  public :: multi_forbar
!
  type forbar
    private
    class(word_iterator), allocatable :: iter
  contains
    procedure :: running => forbar_running
    final :: forbar_destroy
  end type forbar
!
  type multi_forbar
    private
    type(forbar), allocatable :: f(:)
  contains
    procedure :: add     => multi_forbar_add
    procedure :: running => multi_forbar_running
    final :: multi_forbar_destroy
  end type multi_forbar
!
  interface forbar
    module procedure forbar_new
  end interface forbar
!
  interface multi_forbar
    module procedure multi_forbar_new
  end interface multi_forbar
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
  pure function multi_forbar_new() result(res)
    type(multi_forbar) :: res
    ALLOCATE(res%f(0))
  end function multi_forbar_new
!
!| add forbar
  pure subroutine multi_forbar_add(this, iter)
    class(multi_forbar), intent(inout) :: this
    class(word_iterator), intent(in)   :: iter
    type(forbar), allocatable          :: f(:)
    integer                            :: i, n
    n = SIZE(this%f)
    allocate (f(n + 1))
    do concurrent(i=1:n)
      f(i) = this%f(i)
    end do
    f(n + 1) = forbar(iter)
    call MOVE_ALLOC(from=f, to=this%f)
  end subroutine multi_forbar_add
!
  function forbar_running(this) result(res)
    class(forbar), intent(inout) :: this
    logical                      :: res
      call this%iter%next()
      write (*, '(A)', advance='NO') this%iter%var
      res = this%iter%running()
  end function forbar_running
!
  function multi_forbar_running(this, unit) result(res)
    class(multi_forbar), intent(inout) :: this
    integer, intent(in), optional      :: unit
    logical                            :: res
    integer                            :: unit_
    integer                            :: i, ios
    logical                            :: to_tty
!
    unit_ = optarg(unit, STDOUT)
    to_tty = isatty(unit_)
!
    do i = 1, SIZE(this%f)
      call this%f(i)%iter%next()
      if (this%f(i)%iter%to_tty() .and. to_tty .or. &
     &   .not. (this%f(i)%iter%to_tty() .or. to_tty) &
     &   ) then
        write (unit_, '(A)', advance='NO', iostat=ios) this%f(i)%iter%var
      endif
    end do
!
    res = .true.
    do i = 1, SIZE(this%f)
      res = res .and. this%f(i)%iter%running()
    end do
!
  end function multi_forbar_running
!
!| destractor
  pure elemental subroutine forbar_destroy(this)
    type(forbar), intent(inout) :: this
    if (ALLOCATED(this%iter)) deallocate (this%iter)
  end subroutine forbar_destroy
!
!| destractor
  pure elemental subroutine multi_forbar_destroy(this)
    type(multi_forbar), intent(inout) :: this
    if (ALLOCATED(this%f)) deallocate (this%f)
  end subroutine multi_forbar_destroy
!
end module mod_forbar

