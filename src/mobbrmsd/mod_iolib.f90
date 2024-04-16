module mod_iolib
  use ISO_C_BINDING
  use ISO_FORTRAN_ENV, only: INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
  implicit none
  private
  public :: mod_iolib_NEWLINE
  public :: mod_iolib_NULLCHR
  public :: mod_iolib_CAREDGE_RETURN
  public :: mod_iolib_ESCAPE
  public :: mod_iolib_CLEAR_LINE_AFTER
  public :: mod_iolib_CLEAR_LINE_BEFORE
  public :: mod_iolib_CLEAR_LINE
  public :: mod_iolib_CLEAR_SCREEN_AFTER
  public :: mod_iolib_CLEAR_SCREEN_BEFORE
  public :: mod_iolib_CLEAR_SCREEN
  public :: mod_iolib_FS_RESET
  public :: mod_iolib_FS_BOLD
  public :: mod_iolib_FS_WEAK
  public :: mod_iolib_FS_UNDER_LINE
  public :: mod_iolib_FS_INVERT
  public :: mod_iolib_FS_CROSSED_OUT
  public :: mod_iolib_FC_BLACK
  public :: mod_iolib_FC_RED
  public :: mod_iolib_FC_GREEN
  public :: mod_iolib_FC_YELLOW
  public :: mod_iolib_FC_BLUE
  public :: mod_iolib_FC_MAGENTA
  public :: mod_iolib_FC_CYAN
  public :: mod_iolib_FC_WHITE
  public :: isatty
  public :: decorate
  public :: decorator
!&<
!| Newline
  character(*), parameter :: mod_iolib_NEWLINE             = NEW_LINE(' ')
!| Null character
  character(*), parameter :: mod_iolib_NULLCHR             = ACHAR(z'00')
!| Caredge return
  character(*), parameter :: mod_iolib_CAREDGE_RETURN      = ACHAR(z'0d')
!| Escape sequence
  character(*), parameter :: mod_iolib_ESCAPE              = ACHAR(z'1b')
!
!| Clear line after cursor
  character(*), parameter :: mod_iolib_CLEAR_LINE_AFTER    = mod_iolib_ESCAPE//'[K'
!| Clear line before cursor
  character(*), parameter :: mod_iolib_CLEAR_LINE_BEFORE   = mod_iolib_ESCAPE//'[1K'
!| Clear line
  character(*), parameter :: mod_iolib_CLEAR_LINE          = mod_iolib_ESCAPE//'[2K'
!
!| Clear screen after cursor
  character(*), parameter :: mod_iolib_CLEAR_SCREEN_AFTER  = mod_iolib_ESCAPE//'[J'
!| Clear screen before cursor
  character(*), parameter :: mod_iolib_CLEAR_SCREEN_BEFORE = mod_iolib_ESCAPE//'[1J'
!| Clear screen
  character(*), parameter :: mod_iolib_CLEAR_SCREEN        = mod_iolib_ESCAPE//'[2J'
!
!| Font style reset
  character(*), parameter :: mod_iolib_FS_RESET            = mod_iolib_ESCAPE//'[0m'
!
!| Font style bold
  character(*), parameter :: mod_iolib_FS_BOLD             = mod_iolib_ESCAPE//'[1m'
!| Font style weak
  character(*), parameter :: mod_iolib_FS_WEAK             = mod_iolib_ESCAPE//'[2m'
!| Font style underline
  character(*), parameter :: mod_iolib_FS_UNDER_LINE       = mod_iolib_ESCAPE//'[4m'
!| Font style inbert style
  character(*), parameter :: mod_iolib_FS_INVERT           = mod_iolib_ESCAPE//'[7m'
!| Font style corssed out
  character(*), parameter :: mod_iolib_FS_CROSSED_OUT      = mod_iolib_ESCAPE//'[9m'
!
!| Font color black
  character(*), parameter :: mod_iolib_FC_BLACK            = mod_iolib_ESCAPE//'[30m'
!| Font color red
  character(*), parameter :: mod_iolib_FC_RED              = mod_iolib_ESCAPE//'[31m'
!| Font color green
  character(*), parameter :: mod_iolib_FC_GREEN            = mod_iolib_ESCAPE//'[32m'
!| Font color yellow
  character(*), parameter :: mod_iolib_FC_YELLOW           = mod_iolib_ESCAPE//'[33m'
!| Font color blue
  character(*), parameter :: mod_iolib_FC_BLUE             = mod_iolib_ESCAPE//'[34m'
!| Font color maggenta
  character(*), parameter :: mod_iolib_FC_MAGENTA          = mod_iolib_ESCAPE//'[35m'
!| Font color cyan
  character(*), parameter :: mod_iolib_FC_CYAN             = mod_iolib_ESCAPE//'[36m'
!| Font color white
  character(*), parameter :: mod_iolib_FC_WHITE            = mod_iolib_ESCAPE//'[37m'
!&>
  interface
    function isatty_stdin() bind(c)
      use ISO_C_BINDING
      integer(C_INT) :: isatty_stdin
    end function isatty_stdin
!
    function isatty_stdout() bind(c)
      use ISO_C_BINDING
      integer(C_INT) :: isatty_stdout
    end function isatty_stdout
!
    function isatty_stderr() bind(c)
      use ISO_C_BINDING
      integer(C_INT) :: isatty_stderr
    end function isatty_stderr
  end interface
!
  interface optarg
    module procedure :: optarg_l, optarg_c
  end interface optarg
!
contains
!
  !| Return true if unit is tty
  function isatty(unit) result(res)
    integer, intent(in), optional :: unit
    !| Input unit
    logical                       :: res
    integer                       :: cres
    if (.not. PRESENT(unit)) then
      res = isatty_stdout() > 0
      return
    end if
    select case (unit)
    case (INPUT_UNIT)
      cres = isatty_stdin()
    case (OUTPUT_UNIT)
      cres = isatty_stdout()
    case (ERROR_UNIT)
      cres = isatty_stderr()
    case default
      cres = 0
    end select
    res = cres > 0
  end function isatty
!
  !| Return decorator string
  pure function decorate(s, color, style, carret, clear, trimed) result(res)
    character(*), intent(in)        :: s
    !| decorated string <br>
    character, intent(in), optional :: color
    !| Color specifier <br>
    !! 'K' :: Black <br>
    !! 'R' :: Red <br>
    !! 'G' :: Green <br>
    !! 'B' :: Brue <br>
    !! 'M' :: Magenta <br>
    !! 'C' :: Cyan <br>
    !! 'Y' :: Yellow <br>
    !! 'W' :: White <br>
    character, intent(in), optional :: style
    !| Style string <br>
    !! 'B' :: Bold <br>
    !! 'W' :: Weak <br>
    !! 'U' :: Under line <br>
    !! 'I' :: Invert <br>
    !! 'C' :: Crossed out <br>
    logical, intent(in), optional   :: carret
    !| If true, caredge return <br>
    character, intent(in), optional :: clear
    !| Clear options (before carret)<br>
    !! 'L' :: Clear line <br>
    !! 'a' :: Clear line after cursor <br>
    !! 'b' :: Clear line before cursor <br>
    !! 'S' :: Clear screen <br>
    !! 'A' :: Clear screen after cursor <br>
    !! 'B' :: Clear screen before cursor <br>
    logical, intent(in), optional   :: trimed
    !| If true, trim s <br>
    character(:), allocatable       :: dec, res
    !| decorator string
    dec = decorator(color, style, carret, clear)
    if(dec=='')then
      if(optarg(trimed, .false.))then
        res = s
      else
        res = TRIM(s)
      endif
    else
      if(optarg(trimed, .false.))then
        res = dec//s//mod_iolib_FS_RESET
      else
        res = dec//TRIM(s)//mod_iolib_FS_RESET
      endif
    endif
  end function decorate
!
  pure function decorator(color, style, carret, clear) result(res)
    character, intent(in), optional :: color
    !| Color specifier <br>
    !! 'K' :: Black <br>
    !! 'R' :: Red <br>
    !! 'G' :: Green <br>
    !! 'B' :: Brue <br>
    !! 'M' :: Magenta <br>
    !! 'C' :: Cyan <br>
    !! 'Y' :: Yellow <br>
    !! 'W' :: White <br>
    character, intent(in), optional :: style
    !| Style string <br>
    !! 'B' :: Bold <br>
    !! 'W' :: Weak <br>
    !! 'U' :: Under line <br>
    !! 'I' :: Invert <br>
    !! 'C' :: Crossed out <br>
    logical, intent(in), optional   :: carret
    !| If true, caredge return <br>
    character, intent(in), optional :: clear
    !| Clear options (before carret)<br>
    !! 'L' :: Clear line <br>
    !! 'a' :: Clear line after cursor <br>
    !! 'b' :: Clear line before cursor <br>
    !! 'S' :: Clear screen <br>
    !! 'A' :: Clear screen after cursor <br>
    !! 'B' :: Clear screen before cursor <br>
    character(:), allocatable       :: res
    !| decorator string
!&<
    select case (optarg(color, ' '))
    case ('k', 'K'); call decorate_(mod_iolib_FC_BLACK,   style, carret, clear, res)
    case ('r', 'R'); call decorate_(mod_iolib_FC_RED,     style, carret, clear, res)
    case ('g', 'G'); call decorate_(mod_iolib_FC_GREEN,   style, carret, clear, res)
    case ('b', 'B'); call decorate_(mod_iolib_FC_BLUE,    style, carret, clear, res)
    case ('m', 'M'); call decorate_(mod_iolib_FC_MAGENTA, style, carret, clear, res)
    case ('c', 'C'); call decorate_(mod_iolib_FC_CYAN,    style, carret, clear, res)
    case ('y', 'Y'); call decorate_(mod_iolib_FC_YELLOW,  style, carret, clear, res)
    case ('w', 'W'); call decorate_(mod_iolib_FC_WHITE,   style, carret, clear, res)
    case default;    call decorate_('',                   style, carret, clear, res)
    end select
!&>
  end function decorator
!
  pure subroutine decorate_(color, style, carret, clear, res)
    character(*), intent(in)                 :: color
    character(*), intent(in), optional       :: style
    logical, intent(in), optional            :: carret
    character, intent(in), optional          :: clear
    character(:), allocatable, intent(inout) :: res
!&<
    select case (optarg(style, ' '))
    case ('b', 'B'); call decorate__(color, mod_iolib_FS_BOLD,        carret, clear, res)
    case ('w', 'W'); call decorate__(color, mod_iolib_FS_WEAK,        carret, clear, res)
    case ('u', 'U'); call decorate__(color, mod_iolib_FS_UNDER_LINE,  carret, clear, res)
    case ('i', 'I'); call decorate__(color, mod_iolib_FS_INVERT,      carret, clear, res)
    case ('c', 'C'); call decorate__(color, mod_iolib_FS_CROSSED_OUT, carret, clear, res)
    case default;    call decorate__(color, '',                       carret, clear, res)
    end select
!&>
  end subroutine decorate_
!
  pure subroutine decorate__(color, style, carret, clear, res)
    character(*), intent(in)                 :: color
    character(*), intent(in)                 :: style
    logical, intent(in), optional            :: carret
    character, intent(in), optional          :: clear
    character(:), allocatable, intent(inout) :: res
!
    if (optarg(carret, .false.)) then
      call decorate___(color, style, mod_iolib_CAREDGE_RETURN, clear, res)
    else
      call decorate___(color, style, '', clear, res)
    end if
!
  end subroutine decorate__
!
  pure subroutine decorate___(color, style, carret, clear, res)
    character(*), intent(in)                 :: color
    character(*), intent(in)                 :: style
    character(*), intent(in)                 :: carret
    character, intent(in), optional          :: clear
    character(:), allocatable, intent(inout) :: res
!&<
    select case (optarg(clear, ' '))
    case ('L', 'l'); call decorate____(color, style, carret, mod_iolib_CLEAR_LINE,          res)
    case ('a');      call decorate____(color, style, carret, mod_iolib_CLEAR_LINE_AFTER,    res)
    case ('b');      call decorate____(color, style, carret, mod_iolib_CLEAR_LINE_BEFORE,   res)
    case ('S', 's'); call decorate____(color, style, carret, mod_iolib_CLEAR_SCREEN,        res)
    case ('A');      call decorate____(color, style, carret, mod_iolib_CLEAR_SCREEN_AFTER,  res)
    case ('B');      call decorate____(color, style, carret, mod_iolib_CLEAR_SCREEN_BEFORE, res)
    case default;    call decorate____(color, style, carret, '',                            res)
    end select
!&>
  end subroutine decorate___
!
  pure subroutine decorate____(color, style, carret, clear, res)
    character(*), intent(in)                 :: color
    character(*), intent(in)                 :: style
    character(*), intent(in)                 :: carret
    character(*), intent(in)                 :: clear
    character(:), allocatable, intent(inout) :: res
!
    res = clear//carret//color//style
!
  end subroutine decorate____
!
  pure function optarg_l(l, def)
    logical, intent(in), optional :: l
    logical, intent(in)           :: def
    logical                       :: optarg_l
    if (PRESENT(l)) then
      optarg_l = l
    else
      optarg_l = def
    end if
  end function optarg_l
!
  pure function optarg_c(l, def)
    character, intent(in), optional :: l
    character, intent(in)           :: def
    character                       :: optarg_c
    if (PRESENT(l)) then
      optarg_c = l
    else
      optarg_c = def
    end if
  end function optarg_c
!
end module mod_iolib
