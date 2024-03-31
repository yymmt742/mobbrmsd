!| Module for molecular coodinate block indicator.<br>
!  Coordinates must be stored in the following format.<br>
!    X(d,m,n)<br>
!    - d :: spatial dimension.<br>
!    - m :: number of atom in a molecule.<br>
!    - n :: number of molecule.
module mod_mol_block
  use blas_lapack_interface, only : D
  use mod_params, only: IK, RK
  use mod_group_permutation
  implicit none
  private
  public :: mol_block
  public :: mol_block_each_size
  public :: mol_block_total_size
  public :: mol_block_natm
  public :: mol_block_nmol
  public :: mol_block_napm
  public :: mol_block_nsym
  public :: mol_block_swap
  public :: mol_block_inverse_swap
!
  integer(IK), parameter  :: m = 1
  integer(IK), parameter  :: n = 2
  integer(IK), parameter  :: g = 3
!
!| molecular block
!  This is mainly used for passing during initialization.
  type mol_block
    sequence
    !| q :: work array.
    integer(IK), allocatable :: q(:)
  end type mol_block
!
  interface mol_block
    module procedure mol_block_new
  end interface mol_block
!
contains
!
!| Constructer
  pure function mol_block_new(m, n, sym) result(res)
    integer(IK), intent(in) :: m
    !! number of molecules
    integer(IK), intent(in) :: n
    !! number of atoms per molecule
    integer(IK), intent(in), optional :: sym(:,:)
    !! symmetric codomains, [[a1,a2,...,am],[b1,b2,...,bm],...].
    type(mol_block)         :: res
    type(group_permutation) :: s
!
    s = group_permutation(sym)
    res%q = [MAX(m, 0), MAX(n, 1), s%q]
!
  end function mol_block_new
!
!| number of atoms
  pure function mol_block_natm(q) result(res)
    integer(IK), intent(in)     :: q(*)
    !! mol_block
    integer(IK)                 :: res
    res = q(m) * q(n)
  end function mol_block_natm
!
!| number of molecules
  pure function mol_block_nmol(q) result(res)
    integer(IK), intent(in)     :: q(*)
    !! mol_block
    integer(IK)                 :: res
    res = q(n)
  end function mol_block_nmol
!
!| number of atoms per molecule
  pure function mol_block_napm(q) result(res)
    integer(IK), intent(in)     :: q(*)
    !! mol_block
    integer(IK)                 :: res
    res = q(m)
  end function mol_block_napm
!
!| number of atoms per molecule
  pure function mol_block_nsym(q) result(res)
    integer(IK), intent(in)     :: q(*)
    !! mol_block
    integer(IK)                 :: res
    res = group_permutation_nsym(q(g))
  end function mol_block_nsym
!
!| memory blocksize per molecule, defined by d*m.
  pure function mol_block_each_size(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! mol_block
    integer(IK)             :: res
    res = D * q(m)
  end function mol_block_each_size
!
!| memory blocksize, defined by d*m*n.
  pure function mol_block_total_size(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! mol_block
    integer(IK)             :: res
    res = D * q(m) * q(n)
  end function mol_block_total_size
!
!| swap.
  pure subroutine mol_block_swap(q, isym, X)
    integer(IK), intent(in) :: q(*)
    !! mol_block
    integer(IK), intent(in) :: isym
    !! symmetry id
    real(RK), intent(inout) :: X(*)
    !! work array
!
    call group_permutation_swap(q(g), isym, D, X)
!
  end subroutine mol_block_swap
!
!| inverse swap.
  pure subroutine mol_block_inverse_swap(q, isym, X)
    integer(IK), intent(in) :: q(*)
    !! mol_block
    integer(IK), intent(in)     :: isym
    !! symmetry id.
    real(RK), intent(inout)     :: X(*)
    !! work array
!
    call group_permutation_inverse(q(g), isym, D, X)
!
  end subroutine mol_block_inverse_swap
!
end module mod_mol_block

