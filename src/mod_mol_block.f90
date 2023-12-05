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
    integer(IK)                  :: d = 0
    integer(IK)                  :: mg = 0
    integer(IK)                  :: mn = 0
    !  d :: spatial dimension
    type(mol_block), allocatable :: b(:)
    !  mol_blocks
  contains
    procedure         :: n_atom   => mol_block_list_n_atom
    procedure         :: n_spc    => mol_block_list_n_spc
    procedure         :: child    => mol_block_list_child
    procedure         :: invalid  => mol_block_list_invalid
    procedure         :: nspecies => mol_block_list_nspecies
    procedure         :: ispecies => mol_block_list_ispecies
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
  pure function mol_block_list_new(d, s, m, n, f) result(res)
    integer(IK), intent(in) :: d
    !  d :: spatial dimension
    integer(IK), intent(in) :: s
    !  s :: number of species
    integer(IK), intent(in) :: m(s)
    !  m :: list of number of atoms in a molecule, m
    integer(IK), intent(in) :: n(s)
    !  n :: list of number of molecules, n
    integer(IK), intent(in) :: f(s)
    !  f :: list of number of free atoms in a molecule, f
    type(mol_block_list)    :: res
    integer(IK)             :: i, p
    res%d = d
    if (s < 1) then; allocate (res%b(0)); return
    else; allocate (res%b(s))
    end if
    p = 1
    do i = 1, s
      res%b(i) = mol_block(p, m(i), n(i), f(i), n(i))
      p = p + d * n(i) * m(i)
    end do
    res%mg = SUM(res%b%m * res%b%g)
    res%mn = SUM(res%b%m * res%b%n)
  end function mol_block_list_new
!
  pure elemental function mol_block_list_n_atom(this) result(res)
    class(mol_block_list), intent(in) :: this
    integer(IK)                       :: res
    if (ALLOCATED(this%b)) then
      res = SUM(this%b%n * this%b%m)
    else
      res = 0
    end if
  end function mol_block_list_n_atom
!
  pure elemental function mol_block_list_n_spc(this) result(res)
    class(mol_block_list), intent(in) :: this
    integer(IK)                       :: res
    if (ALLOCATED(this%b)) then
      res = SIZE(this%b)
    else
      res = 0
    end if
  end function mol_block_list_n_spc
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
      if (res%b(i)%g == 0) cycle
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
  pure elemental function mol_block_invalid(b) result(res)
    type(mol_block), intent(in) :: b
    logical                     :: res
    res = (b%m < 1) .or. (b%n < 1) .or. (b%m < b%f) .or. (b%n < b%g)
  end function mol_block_invalid
!
  pure elemental subroutine mol_block_list_destroy(this)
    type(mol_block_list), intent(inout) :: this
    if (ALLOCATED(this%b)) deallocate (this%b)
  end subroutine mol_block_list_destroy
!
end module mod_mol_block
