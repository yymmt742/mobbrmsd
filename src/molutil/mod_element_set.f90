module mod_element_set
  use mod_params, only: IK, RK, NL => NEWLINE
  use mod_base
  use mod_element
  implicit none
  private
  public :: element_set
!
  integer(IK), parameter  :: ELESET_CHARLEN  = 128
  character(*), parameter :: DEF_title       = 'default_name'
  character(*), parameter :: DEF_description = 'this is default description.'
!
  include 'element_set_CH.f90'
  include 'element_set_heavy.f90'
  include 'element_set_default.f90'
!
  type, extends(base) :: element_set
    private
    character(ELESET_CHARLEN)  :: t = DEF_title
    character(ELESET_CHARLEN)  :: d = DEF_description
    type(element), allocatable :: e(:)
  contains
    procedure :: nlist               => element_set_nlist
    procedure :: title               => element_set_title
    procedure :: description         => element_set_description
    procedure :: formatted_string    => element_set_formatted_string
    procedure :: clear               => element_set_clear
    procedure :: unfmtread           => element_set_unfmtread
    generic   :: read (unformatted)  => unfmtread
    procedure :: unfmtwrite          => element_set_unfmtwrite
    generic   :: write (unformatted) => unfmtwrite
    final     :: element_set_destroy
  end type element_set
!
  interface element_set
    module procedure element_set_new
  end interface element_set
!
contains
!
  pure function element_set_new(setname, ele, title, description) result(res)
    character(*), intent(in), optional  :: setname
    type(element), intent(in), optional :: ele(:)
    character(*), intent(in), optional  :: title
    character(*), intent(in), optional  :: description
    character(2)                        :: setname_
    type(element_set)                   :: res
!
    setname_ = ''
    if(PRESENT(setname))     setname_ = setname
    if(PRESENT(title))       res%t = title
    if(PRESENT(description)) res%d = description
!
    select case (setname_)
    case ('CH', 'ch')
      res%t = title_CH
      res%d = description_CH
      res%e = element_set_CH
    case ('HE', 'he')
      res%t = title_HEAVY
      res%d = description_HEAVY
      res%e = element_set_HEAVY
    case default
      if(PRESENT(ele))then
        res%e = ele
      else
      res%t = title_default
      res%d = description_default
        res%e = element_set_default
      endif
    end select
!
  end function element_set_new
!
  pure elemental function element_set_nlist(this) result(res)
    class(element_set), intent(in) :: this
    integer(IK)                    :: res
    if (ALLOCATED(this%e)) then
      res = SIZE(this%e)
    else
      res = 0
    end if
  end function element_set_nlist
!
  pure function element_set_title(this) result(res)
    class(element_set), intent(in) :: this
    character(:), allocatable      :: res
    res = TRIM(ADJUSTL(this%t))
  end function element_set_title
!
  pure function element_set_description(this) result(res)
    class(element_set), intent(in) :: this
    character(:), allocatable      :: res
    res = TRIM(ADJUSTL(this%d))
  end function element_set_description
!
  pure function element_set_formatted_string(this) result(res)
    use mod_params, only: NL => NEWLINE
    use mod_element, only: EL => ELEMENT_FORMATTED_CHARLEN, &
                    &      EH => ELEMENT_FORMATTED_HEADER, &
                    &      ef => ELEMENT_FORMATTED_STRING
    class(element_set), intent(in) :: this
    character(:), allocatable   :: res
    integer(IK)                 :: l, i, n, h, m
!
    n = element_set_nlist(this)
    m = LEN_TRIM(this%t) + LEN(NL) + &
      & LEN_TRIM(this%d) + LEN(NL) + 2
    h = LEN(EH)
    l = m + EL * n + MERGE(h, 0, n > 0)
!
    allocate (character(l) :: res)
!
    res(:m) = TRIM(this%t)//NL//'- '//TRIM(this%d)//NL
!
    if (n < 1) return
    res(m + 1:m + h) = EH
    do concurrent(i=1:n)
      block
        integer(IK) :: l, u
        u = i * EL + h + m
        l = u - EL + 1
        res(l:u) = ef(this%e(i))
      end block
    end do
  end function element_set_formatted_string
!
  subroutine element_set_unfmtread(this, unit, iostat, iomsg)
    class(element_set), intent(inout) :: this
    integer(IK), intent(in)        :: unit
    integer(IK), intent(out)       :: iostat
    character(*), intent(inout)    :: iomsg
    integer(IK)                    :: n
!
    call element_set_clear(this)
!
    read (unit, IOSTAT=iostat, IOMSG=iomsg) this%t, this%d, n
!
    if (n < 0) return
    allocate (this%e(n))
    if (n == 1) return
    read (unit, IOSTAT=iostat, IOMSG=iomsg) this%e
!
  end subroutine element_set_unfmtread
!
  subroutine element_set_unfmtwrite(this, unit, iostat, iomsg)
    class(element_set), intent(in) :: this
    integer(IK), intent(in)     :: unit
    integer(IK), intent(out)    :: iostat
    character(*), intent(inout) :: iomsg
    integer(IK)                 :: n
    n = element_set_nlist(this)
    if (n > 0) then
      write (unit, IOSTAT=iostat, IOMSG=iomsg) this%t, this%d, n, this%e
    else
      write (unit, IOSTAT=iostat, IOMSG=iomsg) this%t, this%d, n
    end if
  end subroutine element_set_unfmtwrite
!
  pure elemental subroutine element_set_clear(this)
    class(element_set), intent(inout) :: this
    if (ALLOCATED(this%e)) deallocate (this%e)
    this%t = ''
    this%d = ''
  end subroutine element_set_clear
!
  pure elemental subroutine element_set_destroy(this)
    type(element_set), intent(inout) :: this
    call this%clear()
  end subroutine element_set_destroy
!
end module mod_element_set
