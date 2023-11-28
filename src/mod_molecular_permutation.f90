module mod_molecular_permutation
  use mod_params, only: IK, RK
  use mod_optarg
  implicit none
  private
  public :: molecular_permutation
!
  integer(IK), parameter :: DEF_n = 1
  !! default number of molecule.
!
  type :: molecular_permutation
    private
    integer(IK)              :: n = DEF_n
    integer(IK)              :: g = DEF_n
    integer(IK), allocatable :: p(:)
  contains
    procedure :: clear        => molecular_permutation_clear
    final     :: molecular_permutation_destroy
  end type molecular_permutation
!
  interface molecular_permutation
    module procedure molecular_permutation_new
  end interface molecular_permutation
!
contains
!
!| Constructer
  pure function molecular_permutation_new(n, prm) result(res)
    integer(IK), intent(in), optional :: n
    !! default number of molecule = 1.
    integer(IK), intent(in), optional :: prm(:)
    !! permutation indices
    type(molecular_permutation)       :: res
!
    res%n = MAX(optarg(n, DEF_n), 1)
!
    if (PRESENT(sym)) res%sym = sym
    if (PRESENT(fix)) then
      res%fix = fix
      res%free = free_indices(res%n, fix)
    end if
!
  end function molecular_permutation_new
!
  pure function free_indices(n, fix) result(res)
    integer(IK), intent(in)                  :: n, fix(:)
    integer(IK)                              :: res(n - SIZE(fix))
    integer(IK)                              :: i, j, g
    j = 1
    g = SIZE(res)
    do i = 1, n
      if (ANY(i == fix)) cycle
      res(j) = i
      if (j == g) return
      j = j + 1
    end do
  end function free_indices
!
  pure elemental function molecular_permutation_nfix(this) result(res)
    class(molecular_permutation), intent(in) :: this
    integer(IK)                              :: res
    if(ALLOCATED(this%fix))then
      res = SIZE(this%fix)
    else
      res = 0
    endif
  end function molecular_permutation_nfix
!
  pure elemental function molecular_permutation_nfree(this) result(res)
    class(molecular_permutation), intent(in) :: this
    integer(IK)                              :: res
    if(ALLOCATED(this%free))then
      res = SIZE(this%free)
    else
      res = 0
    endif
  end function molecular_permutation_nfree
!
  pure function molecular_permutation_swap_indices(this) result(res)
    class(molecular_permutation), intent(in) :: this
    integer(IK)                              :: res(this%n)
    integer(IK)                              :: i, j, f
!
    if (ALLOCATED(this%fix)) then
      f = this%nfix()
      do concurrent(i=1:f)
        res(i) = this%fix(i)
      end do
      j = f + 1
      do i = 1, this%n
        if (ANY(i == this%fix)) cycle
        res(j) = i
        if (j == this%n) EXIT
        j = j + 1
      end do
    else
      do concurrent(i=1:this%n)
        res(i) = i
      end do
    end if
!
  end function molecular_permutation_swap_indices
!
  pure elemental subroutine molecular_permutation_clear(this)
    class(molecular_permutation), intent(inout) :: this
    if (ALLOCATED(this%sym)) deallocate (this%sym)
    if (ALLOCATED(this%fix)) deallocate (this%fix)
  end subroutine molecular_permutation_clear
!
  pure elemental subroutine molecular_permutation_destroy(this)
    type(molecular_permutation), intent(inout) :: this
    call this%clear()
  end subroutine molecular_permutation_destroy
!
end module mod_molecular_permutation
