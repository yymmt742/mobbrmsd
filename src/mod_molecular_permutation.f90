module mod_molecular_permutation
  use mod_params, only: IK, RK
  use mod_group_permutation
  use mod_molecular_rotation
  implicit none
  private
  public :: molecular_permutation
!
  type, extends(group_permutation) :: molecular_permutation
!   private
    integer(IK)              :: g = 0
    integer(IK), allocatable :: iper(:)
    integer(IK), allocatable :: isym(:)
  contains
    procedure :: n_mol    => molecular_permutation_n_mol
    procedure :: n_free   => molecular_permutation_n_free
    procedure :: child    => molecular_permutation_child
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
  pure elemental function molecular_permutation_new(n) result(res)
    integer(IK), intent(in)           :: n
    !! number of molecules
    type(molecular_permutation)       :: res
    integer(IK)                       :: i
!
    if (n < 0) then
      res%g = 0
      ALLOCATE(res%iper(0))
      ALLOCATE(res%isym(0))
    else
      res%g = n
      ALLOCATE(res%iper(n))
      do concurrent(i=1:n)
        res%iper(i) = i
      end do
      ALLOCATE(res%isym(n), source=0)
    endif
!
    call group_permutation_init(res, res%iper)
!
  end function molecular_permutation_new
!
  pure elemental function molecular_permutation_n_mol(this) result(res)
    class(molecular_permutation), intent(in) :: this
    integer(IK)                              :: res
    if (ALLOCATED(this%iper)) then; res = SIZE(this%iper)
    else; res = 0
    end if
  end function molecular_permutation_n_mol
!
  pure elemental function molecular_permutation_n_free(this) result(res)
    class(molecular_permutation), intent(in) :: this
    integer(IK)                              :: res
    res = this%g
  end function molecular_permutation_n_free
!
  pure elemental function molecular_permutation_child(this, iper, isym) result(res)
    class(molecular_permutation), intent(in) :: this
    integer(IK), intent(in)                  :: iper, isym
    type(molecular_permutation)              :: res
    integer(IK)                              :: i, n
!
    n = this%n_mol()
!
    allocate (res%iper(n))
    allocate (res%isym(n), source=0)
!
    if (.not. ALLOCATED(this%iper)) then
      res%g = 0
      do concurrent(i=1:n)
        res%iper(i) = i
      enddo
    elseif (this%g < 0 .or. iper < 1 .or. this%g < iper) then
      res%g = 0
      do concurrent(i=1:n)
        res%iper(i) = this%iper(i)
      enddo
    else
      res%g = this%g - 1
      do concurrent(i=1:iper-1)
        res%iper(i) = this%iper(i)
      enddo
      do concurrent(i=iper+1:this%g)
        res%iper(i - 1) = this%iper(i)
      enddo
      res%iper(this%g) = this%iper(iper)
      do concurrent(i=this%g+1:n)
        res%iper(i - 1) = this%iper(i)
      enddo
    end if
!
    if (.not. ALLOCATED(this%iper)) then
    else
      do concurrent(i=1:this%g-1)
        res%isym(i) = this%isym(i)
      enddo
      res%isym(this%g) = isym
      do concurrent(i=this%g+1:n)
        res%isym(i) = this%isym(i)
      enddo
    endif
!
    call group_permutation_init(res, res%iper)
!
  end function molecular_permutation_child
!
  pure subroutine molecular_permutation_mol_swap(this, d, m, rot, X)
    class(molecular_permutation), intent(in) :: this
    integer(IK), intent(in)                  :: d, m
    type(molecular_rotation), intent(in)     :: rot
    real(RK), intent(inout)                  :: X(*)
    integer(IK)                              :: i, dm
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
    this%g = 0
    if (ALLOCATED(this%iper)) deallocate (this%iper)
    if (ALLOCATED(this%isym)) deallocate (this%isym)
    call this%clear()
  end subroutine molecular_permutation_destroy
!
!!!
!
! pure function is_injection(n, g, f) result(res)
!   integer(IK), intent(in)           :: n, g
!   integer(IK), intent(in), optional :: f(g)
!   logical                           :: res
!   integer(IK)                       :: i
!   res = .FALSE.
!   if (.not. PRESENT(f)) return
!   do i = 1, g
!     if (f(i) < 1 .or. n < f(i)) return
!   end do
!   do i = 1, g - 1
!     if (ANY(f(i) == f(i + 1:))) return
!   end do
!   res = .TRUE.
! end function is_injection
!
end module mod_molecular_permutation
