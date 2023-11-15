module mod_molecular_permutation
  use mod_params, only: IK, RK
  use mod_optarg
  use mod_molecule
  implicit none
  private
  public :: molecular_permutation
!
  integer(IK), parameter :: DEF_n = 1
  !! default number of molecule = 1.
  integer(IK), parameter :: DEF_d = 3
  !! default spatial dimension = 3.
!
  type :: molecular_permutation
!
    integer(IK)              :: d = 0
    integer(IK)              :: m = 0
    integer(IK)              :: n = 0
    integer(IK), allocatable :: sym(:), fix(:)
!
  contains
!
    procedure :: nfree        => molecular_permutation_nfree
    procedure :: nfix         => molecular_permutation_nfix
    procedure :: free_indices => molecular_permutation_free_indices
    procedure :: swap_indices => molecular_permutation_swap_indices
    procedure :: sort_matrix  => molecular_permutation_sort_matrix
    procedure :: clear        => molecular_permutation_clear
    final     :: molecular_permutation_destroy
!
  end type molecular_permutation
!
  interface molecular_permutation
    module procedure molecular_permutation_new
  end interface molecular_permutation
!
contains
!
!| Generator
  pure function molecular_permutation_new(mol, n, sym, fix, d) result(res)
    class(molecule), intent(in)       :: mol
    !! molecular template.
    integer(IK), intent(in), optional :: n
    !! default number of molecule = 1.
    integer(IK), intent(in), optional :: d
    !! default number of molecule = 1.
    integer(IK), intent(in), optional :: sym(:)
    !! fixed swap
    integer(IK), intent(in), optional :: fix(:)
    !! fixed swap
    type(molecular_permutation)       :: res
!
    res%n = MAX(optarg(n, DEF_n), 1)
    res%m = mol%natom()
    res%d = MAX(optarg(d, DEF_d), 1)
!
    if(PRESENT(sym)) res%sym = sym
    if(PRESENT(fix)) res%fix = fix
!
  end function molecular_permutation_new
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
    res = this%n - this%nfix()
  end function molecular_permutation_nfree
!
  pure function molecular_permutation_free_indices(this) result(res)
    class(molecular_permutation), intent(in) :: this
    integer(IK)                              :: res(this%nfree())
    integer(IK)                              :: i, j, g
    j = 1
    g = this%nfree()
    if (ALLOCATED(this%fix)) then
      do i = 1, this%n
        if (ANY(i == this%fix)) cycle
        res(j) = i
        if (j == g) return
        j = j + 1
      end do
    else
      do concurrent(i=1:this%n)
        res(i) = i
      end do
    end if
  end function molecular_permutation_free_indices
!
  pure function molecular_permutation_swap_indices(this) result(res)
    class(molecular_permutation), intent(in) :: this
    integer(IK)                              :: res(this%n)
    integer(IK)                              :: i, j, f
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
  end function molecular_permutation_swap_indices
!
  pure subroutine molecular_permutation_sort_matrix(this, mol, source, dest)
    class(molecular_permutation), intent(in) :: this
    class(molecule), intent(in)              :: mol
    real(RK), intent(in)                     :: source(*)
    real(RK), intent(inout)                  :: dest(*)
    integer(IK)                              :: s
!
    if(ALLOCATED(this%sym))then
      call sort_matrix(mol, this%d, mol%natom(), this%n, this%sym, this%swap_indices(), source, dest)
    else
      call sort_matrix(mol, this%d, mol%natom(), this%n, [(1, s=1, this%n)], this%swap_indices(), source, dest)
    endif
!
  contains
!
    pure subroutine sort_matrix(mol, d, m, n, sym, swp, source, dest)
      class(molecule), intent(in) :: mol
      integer(IK), intent(in)     :: d, m, n, sym(n), swp(n)
      real(RK), intent(in)        :: source(d, m, n)
      real(RK), intent(inout)     :: dest(d, m, n)
      integer(IK)                 :: i, j, k
!
      do concurrent(k=1:n)
        block
          integer(IK) :: w(m)
          w = mol%sym_index(sym(swp(k)))
          do concurrent(i=1:d, j=1:m)
            dest(i, j, k) = source(i, w(j), swp(k))
          end do
        end block
      end do
!
    end subroutine sort_matrix
!
  end subroutine molecular_permutation_sort_matrix
!
  pure elemental subroutine molecular_permutation_clear(this)
    class(molecular_permutation), intent(inout) :: this
!
    if (ALLOCATED(this%sym)) deallocate (this%sym)
    if (ALLOCATED(this%fix)) deallocate (this%fix)
!
  end subroutine molecular_permutation_clear
!
  pure elemental subroutine molecular_permutation_destroy(this)
    type(molecular_permutation), intent(inout) :: this
!
    call this%clear()
!
  end subroutine molecular_permutation_destroy
!
end module mod_molecular_permutation
