!| molecular coodinate block indicator
module mod_mol_block
  use mod_params, only: IK, RK, ONE => RONE, FOUR => RFOUR, ZERO => RZERO, RHUGE
  implicit none
  private
  public :: mol_block
  public :: mol_block_list
!
!| molecular block indicator
!  atomic coordinates vector must be stored in the following format.
!    X(d,m,n)
!    - d :: spatial dimension.
!    - m :: number of atom in a molecule.
!    - n :: number of molecule.
!    - where X(d,:f,:g)     :: Free rotatable.
!    -       X(d,f+1:,g+1:) :: Fixed.
  type mol_block
    sequence
    integer(IK) :: p = 1
    !  p :: pointer to memory
    integer(IK) :: s = 1
    !  s :: number of molecular symmetry
    integer(IK) :: m = 1
    !  m :: number of atom in a molecule
    integer(IK) :: n = 1
    !  n :: number of molecule
    integer(IK) :: f = 1
    !  f :: number of free atom in a molecule, must be f<=m.
    integer(IK) :: g = 1
    !  g :: number of free molecule, must be g<=n
  end type mol_block
!
  type mol_block_list
    integer(IK)                           :: d = 0
    integer(IK)                           :: mg = 0
    integer(IK)                           :: mn = 0
    !  d :: spatial dimension
    type(mol_block), allocatable          :: b(:)
    !  mol_blocks
  contains
    procedure         :: child       => mol_block_list_child
    procedure         :: invalid     => mol_block_list_invalid
    procedure         :: natom       => mol_block_list_natom
    procedure         :: nspatial    => mol_block_list_nspatial
    procedure         :: nspecies    => mol_block_list_nspecies
    procedure         :: ispecies    => mol_block_list_ispecies
    procedure         :: ipointer    => mol_block_list_ipointer
    procedure         :: n_res       => mol_block_list_n_res
    procedure         :: res_pointer => mol_block_list_res_pointer
    procedure         :: has_child   => mol_block_list_has_child
    procedure         :: clear       => mol_block_list_clear
    final             :: mol_block_list_destroy
  end type mol_block_list
!
  interface mol_block_list
    module procedure mol_block_list_new
  end interface mol_block_list
!
contains
!
! Constructer
  pure function mol_block_list_new(d, l, b) result(res)
    integer(IK), intent(in)     :: d
    !  d :: spatial dimension
    integer(IK), intent(in)     :: l
    !  s :: number of species
    type(mol_block), intent(in) :: b(l)
    !  molecular block
    type(mol_block_list)        :: res
    integer(IK)                 :: i, p
    res%d = d
    if (l < 1) then
      allocate (res%b(0))
      return
    end if
    allocate (res%b(l))
    do concurrent(i=1:l)
      res%b(i) = b(i)
      res%b(i)%s = MAX(1, res%b(i)%s)
      res%b(i)%m = MAX(0, res%b(i)%m)
      res%b(i)%n = MAX(0, res%b(i)%n)
      res%b(i)%f = MAX(0, MIN(res%b(i)%m, res%b(i)%f))
      res%b(i)%g = MAX(0, MIN(res%b(i)%n, res%b(i)%g))
    end do
    p = 1
    do i = 1, l
      res%b(i)%p = p
      p = p + d * res%b(i)%n * res%b(i)%m
    end do
    res%mg = SUM(res%b%m * res%b%g)
    res%mn = SUM(res%b%m * res%b%n)
  end function mol_block_list_new
!
  pure elemental function mol_block_list_natom(this) result(res)
    class(mol_block_list), intent(in) :: this
    integer(IK)                       :: res
    res = this%mn
  end function mol_block_list_natom
!
  pure elemental function mol_block_list_nspatial(this) result(res)
    class(mol_block_list), intent(in) :: this
    integer(IK)                       :: res
    res = this%d
  end function mol_block_list_nspatial
!
  pure elemental function mol_block_list_child(b) result(res)
    class(mol_block_list), intent(in) :: b
    type(mol_block_list)              :: res
    integer(IK)                       :: i
    res%d = b%d
    res%mg = b%mg
    res%mn = b%mn
    if (.not. ALLOCATED(b%b)) then
      allocate (res%b(0))
      return
    end if
    res%b = b%b
    do i = 1, SIZE(res%b)
      if (res%b(i)%g == 0) cycle
      res%b(i)%g = res%b(i)%g - 1
      return
    end do
  end function mol_block_list_child
!
  pure elemental function mol_block_list_ispecies(b) result(res)
    class(mol_block_list), intent(in) :: b
    integer(IK)                       :: i, res
    if (.not. ALLOCATED(b%b)) then
      res = 0
    else
      do i = 1, SIZE(b%b)
        if (b%b(i)%g == 0) cycle
        res = i; return
      end do
    end if
  end function mol_block_list_ispecies
!
  pure elemental function mol_block_list_ipointer(b, imol) result(res)
    class(mol_block_list), intent(in) :: b
    integer(IK), intent(in)           :: imol
    integer(IK)                       :: i, res
    res = 0
    if (ALLOCATED(b%b)) then
      do i = 1, SIZE(b%b)
        if (b%b(i)%g == 0) cycle
        res = b%b(i)%p
        if (imol < 1 .or. b%b(i)%n < imol) return
        res = res + (imol - 1) * b%d * b%b(i)%m
        return
      end do
    end if
  end function mol_block_list_ipointer
!
  pure function mol_block_list_n_res(b) result(res)
    class(mol_block_list), intent(in) :: b
    integer(IK)                       :: i, s, res(b%nspecies())
    res = 0
    s = b%nspecies()
    do i = 1, s
      res(i) = (b%b(i)%n - b%b(i)%g) * b%b(i)%m
    end do
  end function mol_block_list_n_res
!
  pure elemental function mol_block_list_res_pointer(b, ispc) result(res)
    class(mol_block_list), intent(in) :: b
    integer(IK), intent(in)           :: ispc
    integer(IK)                       :: res
    res = 0
    if (ALLOCATED(b%b)) then
      if (ispc < 1 .or. SIZE(b%b) < ispc) return
      if (b%b(ispc)%n <= b%b(ispc)%g) return
      res = b%b(ispc)%p + b%d * b%b(ispc)%m * (b%b(ispc)%n - b%b(ispc)%g); return
    end if
  end function mol_block_list_res_pointer
!
  pure elemental function mol_block_list_nspecies(this) result(res)
    class(mol_block_list), intent(in) :: this
    integer(IK)                       :: res
    if (.not. ALLOCATED(this%b)) then
      res = 0
    else
      res = SIZE(this%b)
    end if
  end function mol_block_list_nspecies
!
  pure elemental function mol_block_list_invalid(this) result(res)
    class(mol_block_list), intent(in) :: this
    logical                           :: res
    if (this%d < 1 .or. .not. ALLOCATED(this%b)) then
      res = .false.
    else
      res = ANY(mol_block_invalid(this%b))
    end if
  end function mol_block_list_invalid
!
  pure elemental function mol_block_list_has_child(this) result(res)
    class(mol_block_list), intent(in) :: this
    logical                           :: res
    res = .false.
    if (ALLOCATED(this%b)) res = ANY(this%b%g > 0)
  end function mol_block_list_has_child
!
  pure elemental function mol_block_invalid(b) result(res)
    type(mol_block), intent(in) :: b
    logical                     :: res
    res = (b%m < 1) .or. (b%n < 1) .or. (b%m < b%f) .or. (b%n < b%g)
  end function mol_block_invalid
!
  pure elemental subroutine mol_block_list_clear(this)
    class(mol_block_list), intent(inout) :: this
    if (ALLOCATED(this%b)) deallocate (this%b)
  end subroutine mol_block_list_clear
!
  pure elemental subroutine mol_block_list_destroy(this)
    type(mol_block_list), intent(inout) :: this
    if (ALLOCATED(this%b)) deallocate (this%b)
  end subroutine mol_block_list_destroy
!
end module mod_mol_block
