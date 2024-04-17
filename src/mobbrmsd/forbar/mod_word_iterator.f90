module mod_progress_bar
  use mod_optarg
  use mod_params, only: &
 &      IK, &
 &      RK, &
 &      ONE => RONE, &
 &      ZERO => RZERO, &
 &      RHUGE
  use mod_iolib, only: &
 &      STDOUT => mod_iolib_STDOUT, &
 &      NULL_CHAR => mod_iolib_NULLCHR, &
 &      isatty, &
 &      decorate, &
 &      decorator
  implicit none
  private
  public :: progress_bar
!
!| Unicode BLOCK characters
  character(*), parameter :: FULL_BLOCK           = '█'
  character(*), parameter :: SEVEN_EIGHTHS_BLOCK  = '▉'
  character(*), parameter :: THREE_QUARTERS_BLOCK = '▊'
  character(*), parameter :: FIVE_EIGHTHS_BLOCK   = '▋'
  character(*), parameter :: HALF_BLOCK           = '▌'
  character(*), parameter :: ONE_QUARTER_BLOCK    = '▍'
  character(*), parameter :: THREE_EIGHTHS_BLOCK  = '▎'
  character(*), parameter :: ONE_EIGHTH_BLOCK     = '▏'
!
!| Unicode character SHADE BLOCK characters
  character(*), parameter :: LIGHT_SHADE  = '░'
  character(*), parameter :: MEDIUM_SHADE = '▒'
  character(*), parameter :: DARK_SHADE   = '▓'
!
!| Unicode character LIGHT_SHADE
  integer(IK), parameter  :: LEN_FULL_BLOCK           = LEN(HALF_BLOCK)
  integer(IK), parameter  :: LEN_SEVEN_EIGHTHS_BLOCK  = LEN(SEVEN_EIGHTHS_BLOCK)
  integer(IK), parameter  :: LEN_THREE_QUARTERS_BLOCK = LEN(THREE_QUARTERS_BLOCK)
  integer(IK), parameter  :: LEN_FIVE_EIGHTHS_BLOCK   = LEN(FIVE_EIGHTHS_BLOCK)
  integer(IK), parameter  :: LEN_HALF_BLOCK           = LEN(HALF_BLOCK)
  integer(IK), parameter  :: LEN_ONE_QUARTER_BLOCK    = LEN(ONE_QUARTER_BLOCK)
  integer(IK), parameter  :: LEN_THREE_EIGHTHS_BLOCK  = LEN(THREE_EIGHTHS_BLOCK)
  integer(IK), parameter  :: LEN_ONE_EIGHTH_BLOCK     = LEN(ONE_EIGHTH_BLOCK)
  integer(IK), parameter  :: LEN_LIGHT_SHADE          = LEN(LIGHT_SHADE)
  integer(IK), parameter  :: LEN_MEDIUM_SHADE         = LEN(MEDIUM_SHADE)
  integer(IK), parameter  :: LEN_DARK_SHADE           = LEN(DARK_SHADE)
!
  integer(IK), parameter :: LEN_BLOCK = MAX(&
                           &              LEN_FULL_BLOCK, &
                           &              LEN_SEVEN_EIGHTHS_BLOCK, &
                           &              LEN_THREE_QUARTERS_BLOCK, &
                           &              LEN_FIVE_EIGHTHS_BLOCK, &
                           &              LEN_HALF_BLOCK, &
                           &              LEN_ONE_QUARTER_BLOCK, &
                           &              LEN_THREE_EIGHTHS_BLOCK, &
                           &              LEN_ONE_EIGHTH_BLOCK &
                           &            )
!
  character(LEN_BLOCK), parameter :: BLOCKS(8)   = [ &
                                   &  ONE_EIGHTH_BLOCK, &
                                   &  THREE_EIGHTHS_BLOCK, &
                                   &  ONE_QUARTER_BLOCK, &
                                   &  HALF_BLOCK, &
                                   &  FIVE_EIGHTHS_BLOCK, &
                                   &  THREE_QUARTERS_BLOCK, &
                                   &  SEVEN_EIGHTHS_BLOCK, &
                                   &  FULL_BLOCK &
                                   & ]
!
!| NUMBER of ANIMATION frames
  integer(IK), parameter  :: N_ANIMATION_A = 6
!
!| Character length of ANIMATION frame
  integer(IK), parameter  :: L_ANIMATION_A = MAX( &
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
  integer(IK), parameter  :: N_ANIMATION_B = 10
!
!| Character length of ANIMATION frame
  integer(IK), parameter  :: L_ANIMATION_B = MAX( &
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
  type progress_bar
    private
    integer(IK)               :: limit
    integer(IK)               :: pprev
    integer(IK)               :: qprev
    integer(IK)               :: cprev
    integer(IK)               :: counter
    integer(IK)               :: anime_speed
    integer(IK)               :: bar_length
    integer(IK)               :: mid_length
    integer(IK)               :: p1
    integer(IK)               :: p2
    integer(IK)               :: nanim
    character(:), allocatable :: bar_tty
    character(:), allocatable :: layer
    character(:), allocatable :: ifmt
    character(:), allocatable :: cmid
    character(:), allocatable :: deco
    character(:), allocatable :: dres
  contains
    procedure :: update => progress_bar_update
    final :: progress_bar_destroy
  end type progress_bar
!
  interface progress_bar
    module procedure progress_bar_new
  end interface progress_bar
!
contains
!
  pure function num_digit(s) result(res)
    integer(IK), intent(in) :: s ! must be natural number
    integer(IK)             :: res
    res = INT(LOG10(real(MAX(s, 1), RK)), IK) + 2
  end function num_digit
!
!| constructor
  pure function progress_bar_new(limit, title, indent, bar_length) result(res)
    integer(IK), intent(in)           :: limit
    !! mobbrmsd_input
    character(*), intent(in)          :: title
    !! mobbrmsd_input
    integer(IK), intent(in), optional :: indent
    !! indent depth
    integer(IK), intent(in), optional :: bar_length
    !! bar_length
    type(progress_bar)                :: res
!
    integer(IK)                       :: p3
    integer(IK)                       :: indent_
    integer(IK)                       :: digit, digres
    character(:), allocatable         :: header, footer
    character(32)                     :: ifmt
!   character(:), allocatable         :: deco_reset
!   character(:), allocatable         :: deco_end
!
    indent_ = MAX(optarg(indent, 0), 0)
    res%bar_length = MIN(MAX(optarg(bar_length, 60) - indent_, 40), 80)
    res%limit = MAX(limit, 1)
    res%pprev = 0
    res%qprev = 0
    res%cprev = 0
    res%counter = 0
    res%anime_speed = res%limit / 25
!
    res%nanim = N_ANIMATION_B
!
    digit = num_digit(res%limit)
    digres = res%bar_length - digit - LEN(' of') - LEN(' ')
    write (ifmt, '(A,I0,A,I0,A)') '(I', digres, '," of",I', digit, '," ")'
    res%ifmt = TRIM(ifmt)
!
    allocate (character(res%bar_length) :: res%layer)
!
    header = decorator(carret=.true., clear='a')// &
           & REPEAT(' ', indent_)// &
           & decorate( &
           &   title, &
           &   style='Underline' &
           & )// &
           & ' ▍▏'// &
           & decorator(color='K', bgcolor='Green', reset=.true.)
    footer = decorator(reset=.true.)//' ▏▍ '
!
    res%deco = decorator(color='Green', bgcolor='Red', reset=.true.)
    res%dres = decorator(bgcolor='Red', reset=.true.)
!
    res%mid_length = &
        &    LEN(res%deco) + &
        &    LEN_BLOCK + &
        &    LEN(res%dres)
!
    allocate (character(res%mid_length) :: res%cmid)
    call composite_mid(res%deco, BLOCKS(1), res%dres, res%cmid)
!
    res%p1 = LEN(header) + 1
    res%p2 = res%p1 - 1 + res%bar_length + res%mid_length
    p3 = res%p2 + LEN(footer)
!
    allocate (character(p3) :: res%bar_tty)
!
    res%bar_tty(:res%p1-1) = header
    call composite_bar( &
   &       res%bar_length, 0, res%limit, 0, &
   &       res%cmid, res%ifmt, res%layer, &
   &       res%bar_tty(res%p1:res%p2))
    res%bar_tty(res%p2+1:) = footer
!
  end function progress_bar_new
!
  pure subroutine composite_mid(deco, b, dres, res)
    character(*), intent(in)    :: deco, b, dres
    character(*), intent(inout) :: res
    integer(IK)                 :: i, p1, p2
    p1 = 1
    p2 = LEN(deco) + LEN(b) + LEN(dres)
    res(p1:p2) = deco//b//dres
    p1 = p2 + 1
    p2 = LEN(res)
    do concurrent(i=p1:p2)
      res(i:i) = NULL_CHAR
    end do
  end subroutine composite_mid
!
  pure subroutine composite_bar(l, p, limit, prog, cmid, ifmt, layer, res)
    integer(IK), intent(in)     :: l, p, limit, prog
    character(*), intent(in)    :: cmid
    character(*), intent(in)    :: ifmt
    character(*), intent(inout) :: layer
    character(*), intent(inout) :: res
    integer(IK)                 :: i, p1, p2, q1, q2
!
    write (layer, ifmt) prog, limit
!
    p1 = 1
    p2 = p
    do concurrent(i=p1:p2)
      res(i:i) = layer(i:i)
    end do
!
    if (p < l) then
      q1 = p2 + 1
      q2 = p2 + LEN(cmid)
      if (q2 <= LEN(res)) res(q1:q2) = cmid
!
      p1 = p2 + 1
      p2 = l
      q1 = q2 + 1
      q2 = q2 + l - p
      if (q2 <= LEN(res)) res(q1:q2) = layer(p1:p2)
    else
      q2 = p
    end if
!
    q1 = q2 + 1
    q2 = LEN(res)
    do concurrent(i=q1:q2)
      res(i:i) = NULL_CHAR
    end do
!
  end subroutine composite_bar
!
!| update progress_bar
  subroutine progress_bar_update(this, progress)
    class(progress_bar), intent(inout) :: this
    !! this
    integer(IK), intent(in)            :: progress
    !! progress
    integer(IK)                        :: prog_
    integer(IK)                        :: lp, p, q
!
    this%counter = MODULO(this%counter + 1, this%anime_speed)
!
    if (this%counter == 0) then
      this%cprev = MODULO(this%cprev + 1, this%nanim)
    end if
!
    prog_ = MIN(MAX(0, progress), this%limit)
!
    lp = this%bar_length * prog_
    p = lp / this%limit
    q = MODULO(lp, this%limit) * 8 / this%limit + 1
!
    if (this%pprev < p.or.this%qprev < q) then
      this%pprev = p
      this%qprev = q
!
      call composite_mid(this%deco, BLOCKS(q), this%dres, this%cmid)
      call composite_bar(this%bar_length, p, this%limit, prog_, &
     &                   this%cmid, this%ifmt, this%layer, &
     &                   this%bar_tty(this%p1:this%p2))
!
      write (*, '(A)', ADVANCE='NO') this%bar_tty
      FLUSH (6)
!
    end if
!
  end subroutine progress_bar_update
!
!| destractor
  pure elemental subroutine progress_bar_destroy(this)
    type(progress_bar), intent(inout) :: this
    if (ALLOCATED(this%bar_tty)) deallocate (this%bar_tty)
    if (ALLOCATED(this%layer))   deallocate (this%layer)
    if (ALLOCATED(this%ifmt))    deallocate (this%ifmt)
    if (ALLOCATED(this%deco))    deallocate (this%deco)
    if (ALLOCATED(this%dres))    deallocate (this%dres)
  end subroutine progress_bar_destroy
!
end module mod_progress_bar

