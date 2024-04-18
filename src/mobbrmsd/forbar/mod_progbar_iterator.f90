module mod_progbar_iterator
  use mod_optarg
  use mod_iolib
  use mod_word_iterator
  implicit none
  private
  public :: progbar_iterator
!&<
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
!| Unicode character LIGHT_SHADE
  integer, parameter  :: LEN_FULL_BLOCK           = LEN(HALF_BLOCK)
  integer, parameter  :: LEN_SEVEN_EIGHTHS_BLOCK  = LEN(SEVEN_EIGHTHS_BLOCK)
  integer, parameter  :: LEN_THREE_QUARTERS_BLOCK = LEN(THREE_QUARTERS_BLOCK)
  integer, parameter  :: LEN_FIVE_EIGHTHS_BLOCK   = LEN(FIVE_EIGHTHS_BLOCK)
  integer, parameter  :: LEN_HALF_BLOCK           = LEN(HALF_BLOCK)
  integer, parameter  :: LEN_ONE_QUARTER_BLOCK    = LEN(ONE_QUARTER_BLOCK)
  integer, parameter  :: LEN_THREE_EIGHTHS_BLOCK  = LEN(THREE_EIGHTHS_BLOCK)
  integer, parameter  :: LEN_ONE_EIGHTH_BLOCK     = LEN(ONE_EIGHTH_BLOCK)
!
  integer, parameter  :: N_BLOCKS  = 8
!
  integer, parameter  :: LEN_BLOCK = MAX(&
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
!&>
!
  type, extends(word_iterator) :: progbar_iterator
    private
    integer                   :: counter
    integer                   :: n_bar, limit
    integer                   :: mid_p1, mid_p2, mid_len
    integer                   :: bar_p1, bar_p2, bar_len
    character(:), allocatable :: layer
    character(:), allocatable :: cmid
    character(:), allocatable :: ifmt
  contains
    procedure :: next => progbar_iterator_next
    final :: progbar_iterator_destroy
  end type progbar_iterator
!
  interface progbar_iterator
    module procedure progbar_iterator_new
  end interface progbar_iterator
!
contains
!
!| constructor
  pure function progbar_iterator_new(&
 &                limit,&
 &                title, &
 &                indent,&
 &                bar_length,&
 &                maincolor,&
 &                subcolor,&
 &                mainfontcolor,&
 &                subfontcolor&
 &              ) result(res)
    integer, intent(in)                :: limit
    !! limit
    character(*), intent(in), optional :: title
    !! title
    integer, intent(in), optional      :: indent
    !! indent depth
    integer, intent(in), optional      :: bar_length
    !! bar_length
    character(*), intent(in), optional :: maincolor
    !! main color, default ['G'reen]
    character(*), intent(in), optional :: subcolor
    !! sub color, default ['R'reen]
    character(*), intent(in), optional :: mainfontcolor
    !! main font color, default ['K']
    character(*), intent(in), optional :: subfontcolor
    !! sub font color, default ['W']
    character(:), allocatable          :: title_, header, footer, deco, dres
    character                          :: main_, sub_, mainf_, subf_
    character(32)                      :: ifmt
    integer                            :: digit, digres, indent_
    type(progbar_iterator)             :: res
!
    res%counter = 0
    res%limit = MAX(limit, 1)
!
    if (PRESENT(title)) then
      allocate (title_, source=title)
    else
      allocate (character(0) :: title_)
    end if
!
    digit = num_digit(res%limit)
    indent_ = MAX(optarg(indent, 0), 0)
    res%n_bar = MAX( &
              &   MIN( &
              &    optarg(bar_length, 60) - indent_ - LEN(title_), &
              &    80 - indent_ - LEN(title_) &
              &   ), &
              &   digit * 2 + LEN(' of') + LEN(' ') &
              & )
    digres = res%n_bar - digit - LEN(' of') - LEN(' ')
!
    write (ifmt, '(A,I0,A,I0,A)') '(I', digres, '," of",I', digit, '," ")'
    res%ifmt = TRIM(ifmt)
!
    main_  = optarg(maincolor, 'G')
    mainf_ = optarg(mainfontcolor, 'K')
    sub_   = optarg(subcolor, 'R')
    subf_  = optarg(subfontcolor, 'W')
!
    header = REPEAT(' ', indent_)// &
           & decorate( &
           &   title_, &
           &   style='Underline' &
           & )// &
           & ' ▍▏'// &
           & decorator( &
           &   color=mainf_, &
           &   bgcolor=main_,  &
           &   reset=.true.)
!
    footer = decorator(reset=.true.)//' ▏▍'
!
    allocate (character(res%n_bar) :: res%layer)
    write (res%layer, res%ifmt) 0, res%limit
!
    deco = decorator(color=main_, bgcolor=sub_, reset=.true.)
    dres = decorator(color=subf_, bgcolor=sub_, reset=.true.)
!
    res%mid_p1 = LEN(deco) + 1
    res%mid_p2 = res%mid_p1 - 1 + LEN_BLOCK
    res%mid_len = res%mid_p2 + LEN(dres)
!
    allocate (character(res%mid_len) :: res%cmid)
    res%cmid(:res%mid_p1-1) = deco
    res%cmid(res%mid_p1:res%mid_p2) = BLOCKS(1)
    res%cmid(res%mid_p2 + 1:) = dres
!
    res%bar_p1 = LEN(header) + 1
    res%bar_p2 = res%bar_p1 - 1 + res%n_bar + res%mid_len
    res%bar_len = res%bar_p2 + LEN(footer)
    allocate (character(res%bar_len) :: res%var)
    res%var(:res%bar_p1-1) = header
    call composite_bar( &
   &       res%n_bar, res%counter, &
   &       res%cmid, res%layer, &
   &       res%var(res%bar_p1:res%bar_p2))
    res%var(res%bar_p2+1:) = footer
!
  end function progbar_iterator_new
!
!| iterator
  pure subroutine progbar_iterator_next(this)
    class(progbar_iterator), intent(inout) :: this
    integer                                :: lp, p, q
!
    this%counter = this%counter + 1
    if (this%counter > this%limit) this%counter = 0
!
    lp = (this%n_bar + 1) * this%counter
    p = (lp - 1) / this%limit
    if (this%counter == 0) then
      q = 1
    else
      q = MODULO(lp - 1, this%limit) * N_BLOCKS / this%limit + 1
    end if
!
    this%cmid(this%mid_p1:this%mid_p2) = BLOCKS(q)
    write (this%layer, this%ifmt) this%counter, this%limit
    call composite_bar( &
   &       this%n_bar, p, &
   &       this%cmid, this%layer, &
   &       this%var(this%bar_p1:this%bar_p2))
!
  end subroutine progbar_iterator_next
!
  pure subroutine composite_bar(n_bar, bar_counter, cmid, layer, res)
    integer, intent(in)         :: n_bar, bar_counter
    character(*), intent(in)    :: cmid
    character(*), intent(in)    :: layer
    character(*), intent(inout) :: res
    integer                     :: i, p1, p2, q1, q2
!
    p1 = 1
    p2 = bar_counter
    do concurrent(i=p1:p2)
      res(i:i) = layer(i:i)
    end do
!
    q1 = p2 + 1
    q2 = p2 + LEN(cmid)
    if (q2 <= LEN(res)) res(q1:q2) = cmid
!
    p1 = p2 + 1
    p2 = n_bar
    q1 = q2 + 1
    q2 = q2 + n_bar - bar_counter
    if (q2 <= LEN(res)) res(q1:q2) = layer(p1:p2)
!
  end subroutine composite_bar
!
!| destractor
  pure elemental subroutine progbar_iterator_destroy(this)
    type(progbar_iterator), intent(inout) :: this
    if (ALLOCATED(this%var)) deallocate (this%var)
    if (ALLOCATED(this%cmid)) deallocate (this%cmid)
    if (ALLOCATED(this%layer)) deallocate (this%layer)
    if (ALLOCATED(this%ifmt)) deallocate (this%ifmt)
  end subroutine progbar_iterator_destroy
!
  pure function num_digit(s) result(res)
    integer, intent(in) :: s ! must be natural number
    integer             :: res
    res = INT(LOG10(real(MAX(s, 1)))) + 2
  end function num_digit
!
end module mod_progbar_iterator

