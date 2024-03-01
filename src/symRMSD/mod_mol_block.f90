!| Module for molecular coodinate block indicator.<br>
!  Coordinates must be stored in the following format.<br>
!    X(d,m,n)<br>
!    - d :: spatial dimension.<br>
!    - m :: number of atom in a molecule.<br>
!    - n :: number of molecule.
module mod_mol_block
  use mod_params, only: D, IK, RK
  use mod_group_permutation
  implicit none
  private
  public :: mol_block
  public :: mol_block_tuple
  public :: mol_block_set_pointer
  public :: mol_block_pointer
  public :: mol_block_each_size
  public :: mol_block_total_size
  public :: mol_block_nmol
  public :: mol_block_napm
  public :: mol_block_nsym
  public :: mol_block_swap
  public :: mol_block_inverse_swap
!
!| molecular block
  type mol_block
    sequence
    private
    !| p :: pointer to X.
    integer(IK)     :: x = 1
    !| m :: number of atom in a molecule
    integer(IK)     :: m = 1
    !| n :: number of molecule
    integer(IK)     :: n = 1
    !| s :: molecular symmetry
    type(group_permutation) :: s
  end type mol_block
!
  interface mol_block
    module procedure mol_block_new
  end interface mol_block
!
!| A set of molecular block and work arrays. <br>
!  This is mainly used for passing during initialization.
  type mol_block_tuple
    !| m :: number of atom in a molecule
    type (mol_block)         :: b
    !| iw :: work array.
    integer(IK), allocatable :: w(:)
  end type mol_block_tuple
!
  interface mol_block_tuple
    module procedure mol_block_tuple_new
  end interface mol_block_tuple
!
contains
!
  pure function mol_block_tuple_new(m, n, sym) result(res)
    integer(IK), intent(in) :: m
    !! number of molecules
    integer(IK), intent(in) :: n
    !! number of atoms per molecule
    integer(IK), intent(in), optional :: sym(:,:)
    !! symmetric codomains, [[a1,a2,...,am],[b1,b2,...,bm],...].
    type(mol_block_tuple)   :: res
    type(group_permutation_tuple) :: s
!
    res%b = mol_block_new(m, n, sym)
    s = group_permutation_tuple(sym)
    res%w = s%w
!
  end function mol_block_tuple_new
!
! Constructer
  pure function mol_block_new(m, n, sym) result(res)
    integer(IK), intent(in) :: m
    !! number of molecules
    integer(IK), intent(in) :: n
    !! number of atoms per molecule
    integer(IK), intent(in), optional :: sym(:,:)
    !! symmetric codomains, [[a1,a2,...,am],[b1,b2,...,bm],...].
    type(mol_block)         :: res
    res%m = MAX(m, 1)
    res%n = MAX(n, 1)
    res%s = group_permutation(sym)
  end function mol_block_new
!
!| set_pointer.
  pure subroutine mol_block_set_pointer(b, x, s)
    type(mol_block), intent(inout) :: b
    !! mol_block.
    integer(IK), intent(in)        :: x
    !! pointer to coordintate.
    integer(IK), intent(in)        :: s
    !! pointer to mol_symmetry.
!
    b%x = x
    b%s%s = s
!
  end subroutine mol_block_set_pointer
!
!| number of molecules
  pure elemental function mol_block_pointer(b) result(res)
    type(mol_block), intent(in) :: b
    !! mol_block
    integer(IK)                 :: res
    res = b%x
  end function mol_block_pointer
!
!| number of molecules
  pure elemental function mol_block_nmol(b) result(res)
    type(mol_block), intent(in) :: b
    !! mol_block
    integer(IK)                 :: res
    res = b%n
  end function mol_block_nmol
!
!| number of atoms per molecule
  pure elemental function mol_block_napm(b) result(res)
    type(mol_block), intent(in) :: b
    !! mol_block
    integer(IK)                 :: res
    res = b%m
  end function mol_block_napm
!
!| number of atoms per molecule
  pure elemental function mol_block_nsym(b) result(res)
    type(mol_block), intent(in) :: b
    !! mol_block
    integer(IK)                 :: res
    res = b%s%s
  end function mol_block_nsym
!
!| memory blocksize per molecule, defined by d*m.
  pure elemental function mol_block_each_size(b) result(res)
    type(mol_block), intent(in) :: b
    !! mol_block
    integer(IK)                 :: res
    res = D * b%m
  end function mol_block_each_size
!
!| memory blocksize, defined by d*m*n.
  pure elemental function mol_block_total_size(b) result(res)
    type(mol_block), intent(in) :: b
    !! mol_block
    integer(IK)                 :: res
    res = D * b%m * b%n
  end function mol_block_total_size
!
!| swap.
  pure subroutine mol_block_swap(b, isym, ms, X)
    type(mol_block), intent(in) :: b
    !! mol_block.
    integer(IK), intent(in)     :: isym
    !! symmetry id.
    integer(IK), intent(in)     :: ms(*)
    real(RK), intent(inout)     :: X(*)
!
    call group_permutation_swap(b%s, ms, isym, D, X)
!
  end subroutine mol_block_swap
!
!| inverse swap.
  pure subroutine mol_block_inverse_swap(b, isym, ms, X)
    type(mol_block), intent(in) :: b
    !! mol_block.
    integer(IK), intent(in)     :: isym
    !! symmetry id.
    integer(IK), intent(in)     :: ms(*)
    real(RK), intent(inout)     :: X(*)
!
    call group_permutation_inverse(b%s, ms, isym, D, X)
!
  end subroutine mol_block_inverse_swap
!
end module mod_mol_block

