module mod_molecular_permutation
  use mod_params, only: IK, RK
  use mod_group_permutation
  use mod_molecular_rotation
  implicit none
  private
  public :: molecular_permutation
!
  type, extends(group_permutation) :: molecular_permutation
    private
    integer(IK), allocatable :: isym(:)
  contains
    procedure :: mol_swap => molecular_permutation_mol_swap
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
  pure function molecular_permutation_new(n, g, f, p, q) result(res)
    integer(IK), intent(in)           :: n
    !! number of molecules
    integer(IK), intent(in)           :: g
    !! number of free molecules, must be g<=n.
    integer(IK), intent(in), optional :: f(n)
    !! injection n to g
    integer(IK), intent(in), optional :: p(g)
    !! f : injection n to g
    integer(IK), intent(in), optional :: q(g)
    !! f : injection n to g
    integer(IK)                       :: prm(n)
    type(molecular_permutation)       :: res
    integer(IK)                       :: i, j
!
    if (0 < n .and. 0 < g .and. n >= g &
   &    .and. is_injection(n, g, f) &
   &    .and. is_injection(g, g, p)) then
      do concurrent(i=1:g)
        prm(p(f(i))) = i
      end do
      j = g + 1
      do i = 1, n
        if (prm(i) > 0) cycle
        prm(i) = j
        j = j + 1
      end do
    else
      do concurrent(i=1:n)
        prm(i) = i
      end do
    end if
!
    call group_permutation_init(res, prm)
!
    if (PRESENT(q)) then
      res%isym = q
    else
      allocate (res%isym(g), source=0)
    end if
!
  end function molecular_permutation_new
!
  pure subroutine molecular_permutation_mol_swap(this, rot, d, m, X)
    class(molecular_permutation), intent(inout) :: this
    type(molecular_rotation), intent(in)        :: rot
    integer(IK), intent(in)                     :: d, m
    real(RK), intent(inout)                     :: X(*)
    integer(IK)                                 :: i, dm
    dm = d * m
    call group_permutation_swap(this, dm, X)
    do concurrent(i=1:SIZE(this%isym))
      block
        integer(IK) :: j
        j = (i - 1) * dm + 1
        call rot%swap(d, X(j), this%isym(i))
      end block
    end do
  end subroutine molecular_permutation_mol_swap
!
  pure elemental subroutine molecular_permutation_destroy(this)
    type(molecular_permutation), intent(inout) :: this
    if (ALLOCATED(this%isym)) deallocate (this%isym)
    call this%clear()
  end subroutine molecular_permutation_destroy
!
!!!
!
  pure function is_injection(n, g, f) result(res)
    integer(IK), intent(in)           :: n, g
    integer(IK), intent(in), optional :: f(g)
    logical                           :: res
    integer(IK)                       :: i
    res = .FALSE.
    if (.not. PRESENT(f)) return
    do i = 1, g
      if (f(i) < 1 .or. n < f(i)) return
    end do
    do i = 1, g - 1
      if (ANY(f(i) == f(i + 1:))) return
    end do
    res = .TRUE.
  end function is_injection
!
end module mod_molecular_permutation
