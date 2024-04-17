module mod_iolib
  use ISO_C_BINDING
  use ISO_FORTRAN_ENV, only: &
 &   mod_iolib_STDIN => INPUT_UNIT, &
 &   mod_iolib_STDOUT => OUTPUT_UNIT, &
 &   mod_iolib_STDERR => ERROR_UNIT
  use mod_optarg
  implicit none
  private
  public :: mod_iolib_STDIN
  public :: mod_iolib_STDOUT
  public :: mod_iolib_STDERR
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
contains
!
  !| Return true if unit is tty
  function isatty(unit) result(res)
    integer, intent(in), optional :: unit
    !| Input unit
    logical                       :: res
    integer                       :: cres
    if (.not. PRESENT(unit)) then
      cres = isatty_stdout()
    else
      select case (unit)
      case (mod_iolib_STDIN)
        cres = isatty_stdin()
      case (mod_iolib_STDOUT)
        cres = isatty_stdout()
      case (mod_iolib_STDERR)
        cres = isatty_stderr()
      case default
        cres = 0
      end select
    end if
    res = .not. cres == 0
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
  pure function decorator(color, style, carret, clear, clear_screen) result(res)
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
    !! 'Z' :: Clear all <br>
    !! 'A' :: Clear after cursor <br>
    !! 'B' :: Clear before cursor <br>
    logical, intent(in), optional   :: clear_screen
    !| If true, clear screen. <br>
    !! else, clear line (default). <br>
    character(:), allocatable       :: res
    !| decorator string
    associate (BK => mod_iolib_FC_BLACK, &
   &           RD => mod_iolib_FC_RED, &
   &           GR => mod_iolib_FC_GREEN, &
   &           BL => mod_iolib_FC_BLUE, &
   &           MG => mod_iolib_FC_MAGENTA, &
   &           CY => mod_iolib_FC_CYAN, &
   &           YL => mod_iolib_FC_YELLOW, &
   &           WH => mod_iolib_FC_WHITE &
   &     )
      select case (optarg(color, ' '))
      case ('k', 'K'); call deco_(BK, style, carret, clear, clear_screen, res)
      case ('r', 'R'); call deco_(RD, style, carret, clear, clear_screen, res)
      case ('g', 'G'); call deco_(GR, style, carret, clear, clear_screen, res)
      case ('b', 'B'); call deco_(BL, style, carret, clear, clear_screen, res)
      case ('m', 'M'); call deco_(MG, style, carret, clear, clear_screen, res)
      case ('c', 'C'); call deco_(CY, style, carret, clear, clear_screen, res)
      case ('y', 'Y'); call deco_(YL, style, carret, clear, clear_screen, res)
      case ('w', 'W'); call deco_(WH, style, carret, clear, clear_screen, res)
      case default; call deco_('', style, carret, clear, clear_screen, res)
      end select
    end associate
  end function decorator
!
  pure subroutine deco_(color, style, carret, clear, clear_screen, res)
    character(*), intent(in)                 :: color
    character(*), intent(in), optional       :: style
    logical, intent(in), optional            :: carret
    character, intent(in), optional          :: clear
    logical, intent(in), optional            :: clear_screen
    character(:), allocatable, intent(inout) :: res
    associate (BO => mod_iolib_FS_BOLD, &
   &           WE => mod_iolib_FS_WEAK, &
   &           UL => mod_iolib_FS_UNDER_LINE, &
   &           IV => mod_iolib_FS_INVERT, &
   &           CO => mod_iolib_FS_CROSSED_OUT &
   &     )
      select case (optarg(style, ' '))
      case ('b', 'B'); call deco__(color, BO, carret, clear, clear_screen, res)
      case ('w', 'W'); call deco__(color, WE, carret, clear, clear_screen, res)
      case ('u', 'U'); call deco__(color, UL, carret, clear, clear_screen, res)
      case ('i', 'I'); call deco__(color, IV, carret, clear, clear_screen, res)
      case ('c', 'C'); call deco__(color, CO, carret, clear, clear_screen, res)
      case default;    call deco__(color, '', carret, clear, clear_screen, res)
      end select
    end associate
  end subroutine deco_
!
  pure subroutine deco__(color, style, carret, clear, clear_screen, res)
    character(*), intent(in)                 :: color
    character(*), intent(in)                 :: style
    logical, intent(in), optional            :: carret
    character, intent(in), optional          :: clear
    logical, intent(in), optional            :: clear_screen
    character(:), allocatable, intent(inout) :: res
!
    associate (CR => mod_iolib_CAREDGE_RETURN, &
   &           CL => mod_iolib_CLEAR_LINE, &
   &           LA => mod_iolib_CLEAR_LINE_AFTER, &
   &           LB => mod_iolib_CLEAR_LINE_BEFORE, &
   &           CS => mod_iolib_CLEAR_SCREEN, &
   &           SA => mod_iolib_CLEAR_SCREEN_AFTER, &
   &           SB => mod_iolib_CLEAR_SCREEN_BEFORE &
   &     )
      if (optarg(carret, .false.)) then
        if (optarg(clear_screen, .false.)) then
          select case (optarg(clear, ' '))
          case ('Z', 'z'); call deco___(color, style, CR, CL, res)
          case ('A', 'a'); call deco___(color, style, CR, LA, res)
          case ('B', 'b'); call deco___(color, style, CR, LB, res)
          case default;    call deco___(color, style, CR, '', res)
          end select
        else
          select case (optarg(clear, ' '))
          case ('Z', 'z'); call deco___(color, style, CR, CS, res)
          case ('A', 'a'); call deco___(color, style, CR, SA, res)
          case ('B', 'b'); call deco___(color, style, CR, SB, res)
          case default;    call deco___(color, style, CR, '', res)
          end select
        endif
      else
        if (optarg(clear_screen, .false.)) then
          select case (optarg(clear, ' '))
          case ('Z', 'z'); call deco___(color, style, '', CL, res)
          case ('A', 'a'); call deco___(color, style, '', LA, res)
          case ('B', 'b'); call deco___(color, style, '', LB, res)
          case default;    call deco___(color, style, '', '', res)
          end select
        else
          select case (optarg(clear, ' '))
          case ('Z', 'z'); call deco___(color, style, '', CS, res)
          case ('A', 'a'); call deco___(color, style, '', SA, res)
          case ('B', 'b'); call deco___(color, style, '', SB, res)
          case default;    call deco___(color, style, '', '', res)
          end select
        endif
      end if
    end associate
!
  end subroutine deco__
!
  pure subroutine deco___(color, style, carret, clear, res)
    character(*), intent(in)                 :: color
    character(*), intent(in)                 :: style
    character(*), intent(in)                 :: carret
    character(*), intent(in)                 :: clear
    character(:), allocatable, intent(inout) :: res
!
    res = clear//carret//color//style
!
  end subroutine deco___
!
end module mod_iolib
