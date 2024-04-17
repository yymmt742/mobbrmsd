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
!
!| Font color black
  character(*), parameter :: mod_iolib_BG_BLACK            = mod_iolib_ESCAPE//'[40m'
!| Font color red
  character(*), parameter :: mod_iolib_BG_RED              = mod_iolib_ESCAPE//'[41m'
!| Font color green
  character(*), parameter :: mod_iolib_BG_GREEN            = mod_iolib_ESCAPE//'[42m'
!| Font color yellow
  character(*), parameter :: mod_iolib_BG_YELLOW           = mod_iolib_ESCAPE//'[43m'
!| Font color blue
  character(*), parameter :: mod_iolib_BG_BLUE             = mod_iolib_ESCAPE//'[44m'
!| Font color maggenta
  character(*), parameter :: mod_iolib_BG_MAGENTA          = mod_iolib_ESCAPE//'[45m'
!| Font color cyan
  character(*), parameter :: mod_iolib_BG_CYAN             = mod_iolib_ESCAPE//'[46m'
!| Font color white
  character(*), parameter :: mod_iolib_BG_WHITE            = mod_iolib_ESCAPE//'[47m'
!&>
!
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
  pure function decorate( &
 &                s, &
 &                color, &
 &                bgcolor, &
 &                style, &
 &                carret, &
 &                reset, &
 &                clear, &
 &                clear_screen, &
 &                trimed &
 &              ) result(res)
    character(*), intent(in)        :: s
    !| string for decorated <br>
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
    character, intent(in), optional :: bgcolor
    !| Backgoround color specifier <br>
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
    logical, intent(in), optional   :: reset
    !| If true, reset previous font style <br>
    character, intent(in), optional :: clear
    !| Clear options (before carret)<br>
    !! 'Z' :: Clear all <br>
    !! 'A' :: Clear after cursor <br>
    !! 'B' :: Clear before cursor <br>
    logical, intent(in), optional   :: clear_screen
    !| If true, clear screen. <br>
    !! else, clear line (default). <br>
    logical, intent(in), optional   :: trimed
    !| If true, trim s <br>
    character(:), allocatable       :: dec, res
    !| decorator string
    dec = decorator( &
   &        color, &
   &        bgcolor, &
   &        style, &
   &        carret, &
   &        reset, &
   &        clear, &
   &        clear_screen)
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
  pure function decorator( &
 &                color, &
 &                bgcolor, &
 &                style, &
 &                carret, &
 &                reset, &
 &                clear, &
 &                clear_screen &
 &              ) result(res)
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
    character, intent(in), optional :: bgcolor
    !| Backgoround color specifier <br>
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
    logical, intent(in), optional   :: reset
    !| If true, reset previous font style <br>
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
      case ('k', 'K'); call deco_(BK, bgcolor, style, carret, reset, clear, clear_screen, res)
      case ('r', 'R'); call deco_(RD, bgcolor, style, carret, reset, clear, clear_screen, res)
      case ('g', 'G'); call deco_(GR, bgcolor, style, carret, reset, clear, clear_screen, res)
      case ('b', 'B'); call deco_(BL, bgcolor, style, carret, reset, clear, clear_screen, res)
      case ('m', 'M'); call deco_(MG, bgcolor, style, carret, reset, clear, clear_screen, res)
      case ('c', 'C'); call deco_(CY, bgcolor, style, carret, reset, clear, clear_screen, res)
      case ('y', 'Y'); call deco_(YL, bgcolor, style, carret, reset, clear, clear_screen, res)
      case ('w', 'W'); call deco_(WH, bgcolor, style, carret, reset, clear, clear_screen, res)
      case default; call deco_('', bgcolor, style, carret, reset, clear, clear_screen, res)
      end select
    end associate
  end function decorator
!
  pure subroutine deco_(color, bgcolor, style, carret, reset, clear, clear_screen, res)
    character(*), intent(in)                 :: color
    character(*), intent(in), optional       :: bgcolor
    character(*), intent(in), optional       :: style
    logical, intent(in), optional            :: carret
    logical, intent(in), optional            :: reset
    character, intent(in), optional          :: clear
    logical, intent(in), optional            :: clear_screen
    character(:), allocatable, intent(inout) :: res
    associate (BK => mod_iolib_BG_BLACK, &
   &           RD => mod_iolib_BG_RED, &
   &           GR => mod_iolib_BG_GREEN, &
   &           BL => mod_iolib_BG_BLUE, &
   &           MG => mod_iolib_BG_MAGENTA, &
   &           CY => mod_iolib_BG_CYAN, &
   &           YL => mod_iolib_BG_YELLOW, &
   &           WH => mod_iolib_BG_WHITE &
   &     )
      select case (optarg(bgcolor, ' '))
      case ('k', 'K'); call deco__(color, BK, style, carret, reset, clear, clear_screen, res)
      case ('r', 'R'); call deco__(color, RD, style, carret, reset, clear, clear_screen, res)
      case ('g', 'G'); call deco__(color, GR, style, carret, reset, clear, clear_screen, res)
      case ('b', 'B'); call deco__(color, BL, style, carret, reset, clear, clear_screen, res)
      case ('m', 'M'); call deco__(color, MG, style, carret, reset, clear, clear_screen, res)
      case ('c', 'C'); call deco__(color, CY, style, carret, reset, clear, clear_screen, res)
      case ('y', 'Y'); call deco__(color, YL, style, carret, reset, clear, clear_screen, res)
      case ('w', 'W'); call deco__(color, WH, style, carret, reset, clear, clear_screen, res)
      case default; call deco__(color, '', style, carret, reset, clear, clear_screen, res)
      end select
    end associate
  end subroutine deco_
!
  pure subroutine deco__(color, bgcolor, style, carret, reset, clear, clear_screen, res)
    character(*), intent(in)                 :: color
    character(*), intent(in)                 :: bgcolor
    character(*), intent(in), optional       :: style
    logical, intent(in), optional            :: carret
    logical, intent(in), optional            :: reset
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
      case ('b', 'B'); call deco___(color, bgcolor, BO, carret, reset, clear, clear_screen, res)
      case ('w', 'W'); call deco___(color, bgcolor, WE, carret, reset, clear, clear_screen, res)
      case ('u', 'U'); call deco___(color, bgcolor, UL, carret, reset, clear, clear_screen, res)
      case ('i', 'I'); call deco___(color, bgcolor, IV, carret, reset, clear, clear_screen, res)
      case ('c', 'C'); call deco___(color, bgcolor, CO, carret, reset, clear, clear_screen, res)
      case default;    call deco___(color, bgcolor, '', carret, reset, clear, clear_screen, res)
      end select
    end associate
  end subroutine deco__
!
  pure subroutine deco___(color, bgcolor, style, carret, reset, clear, clear_screen, res)
    character(*), intent(in)                 :: color
    character(*), intent(in)                 :: bgcolor
    character(*), intent(in)                 :: style
    logical, intent(in), optional            :: carret
    logical, intent(in), optional            :: reset
    character, intent(in), optional          :: clear
    logical, intent(in), optional            :: clear_screen
    character(:), allocatable, intent(inout) :: res
!
    associate (RE => mod_iolib_FS_RESET, &
   &           CR => mod_iolib_CAREDGE_RETURN, &
   &           CL => mod_iolib_CLEAR_LINE, &
   &           LA => mod_iolib_CLEAR_LINE_AFTER, &
   &           LB => mod_iolib_CLEAR_LINE_BEFORE, &
   &           CS => mod_iolib_CLEAR_SCREEN, &
   &           SA => mod_iolib_CLEAR_SCREEN_AFTER, &
   &           SB => mod_iolib_CLEAR_SCREEN_BEFORE &
   &     )
      if (optarg(reset, .false.)) then
        if (optarg(carret, .false.)) then
          if (optarg(clear_screen, .false.)) then
            select case (optarg(clear, ' '))
            case ('Z', 'z'); call deco____(color, bgcolor, style, CR, RE, CL, res)
            case ('A', 'a'); call deco____(color, bgcolor, style, CR, RE, LA, res)
            case ('B', 'b'); call deco____(color, bgcolor, style, CR, RE, LB, res)
            case default;    call deco____(color, bgcolor, style, CR, RE, '', res)
            end select
          else
            select case (optarg(clear, ' '))
            case ('Z', 'z'); call deco____(color, bgcolor, style, CR, RE, CS, res)
            case ('A', 'a'); call deco____(color, bgcolor, style, CR, RE, SA, res)
            case ('B', 'b'); call deco____(color, bgcolor, style, CR, RE, SB, res)
            case default;    call deco____(color, bgcolor, style, CR, RE, '', res)
            end select
          endif
        else
          if (optarg(clear_screen, .false.)) then
            select case (optarg(clear, ' '))
            case ('Z', 'z'); call deco____(color, bgcolor, style, '', RE, CL, res)
            case ('A', 'a'); call deco____(color, bgcolor, style, '', RE, LA, res)
            case ('B', 'b'); call deco____(color, bgcolor, style, '', RE, LB, res)
            case default;    call deco____(color, bgcolor, style, '', RE, '', res)
            end select
          else
            select case (optarg(clear, ' '))
            case ('Z', 'z'); call deco____(color, bgcolor, style, '', RE, CS, res)
            case ('A', 'a'); call deco____(color, bgcolor, style, '', RE, SA, res)
            case ('B', 'b'); call deco____(color, bgcolor, style, '', RE, SB, res)
            case default;    call deco____(color, bgcolor, style, '', RE, '', res)
            end select
          endif
        end if
      else
        if (optarg(carret, .false.)) then
          if (optarg(clear_screen, .false.)) then
            select case (optarg(clear, ' '))
            case ('Z', 'z'); call deco____(color, bgcolor, style, CR, '', CL, res)
            case ('A', 'a'); call deco____(color, bgcolor, style, CR, '', LA, res)
            case ('B', 'b'); call deco____(color, bgcolor, style, CR, '', LB, res)
            case default;    call deco____(color, bgcolor, style, CR, '', '', res)
            end select
          else
            select case (optarg(clear, ' '))
            case ('Z', 'z'); call deco____(color, bgcolor, style, CR, '', CS, res)
            case ('A', 'a'); call deco____(color, bgcolor, style, CR, '', SA, res)
            case ('B', 'b'); call deco____(color, bgcolor, style, CR, '', SB, res)
            case default;    call deco____(color, bgcolor, style, CR, '', '', res)
            end select
          endif
        else
          if (optarg(clear_screen, .false.)) then
            select case (optarg(clear, ' '))
            case ('Z', 'z'); call deco____(color, bgcolor, style, '', '', CL, res)
            case ('A', 'a'); call deco____(color, bgcolor, style, '', '', LA, res)
            case ('B', 'b'); call deco____(color, bgcolor, style, '', '', LB, res)
            case default;    call deco____(color, bgcolor, style, '', '', '', res)
            end select
          else
            select case (optarg(clear, ' '))
            case ('Z', 'z'); call deco____(color, bgcolor, style, '', '', CS, res)
            case ('A', 'a'); call deco____(color, bgcolor, style, '', '', SA, res)
            case ('B', 'b'); call deco____(color, bgcolor, style, '', '', SB, res)
            case default;    call deco____(color, bgcolor, style, '', '', '', res)
            end select
          endif
        end if
      end if
    end associate
!
  end subroutine deco___
!
  pure subroutine deco____(color, bgcolor, style, carret, reset, clear, res)
    character(*), intent(in)                 :: color
    character(*), intent(in)                 :: bgcolor
    character(*), intent(in)                 :: style
    character(*), intent(in)                 :: carret
    character(*), intent(in)                 :: reset
    character(*), intent(in)                 :: clear
    character(:), allocatable, intent(inout) :: res
!
    res = reset//carret//clear//color//bgcolor//style
!
  end subroutine deco____
!
end module mod_iolib
