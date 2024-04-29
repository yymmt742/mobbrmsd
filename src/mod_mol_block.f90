!| Module for molecular coodinate block indicator.<br>
!  Coordinates must be stored in the 3-D array, X(d, n, M).<br>
!   <br>
!    - \(d\) :: spatial dimension.<br>
!    - \(n\) :: number of atom in a molecule.<br>
!    - \(M\) :: number of molecule.<br>
module mod_mol_block
  use mod_dimspec_functions, only: D
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
  integer(IK), parameter  :: pn = 1
  integer(IK), parameter  :: pm = 2
  integer(IK), parameter  :: pg = 3
!
!| molecular block
!  This is mainly used for passing during initialization.
  type mol_block
    sequence
    !| header
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
  pure function mol_block_new(n, M, sym) result(res)
    integer(IK), intent(in) :: n
    !! number of atoms per molecule
    integer(IK), intent(in) :: M
    !! number of molecules
    integer(IK), intent(in), optional :: sym(:, :)
    !! symmetric codomains, sym(n, S)
    !! \( [[\nu_1^{(1)},\dots,\nu_n^{(1)}],[\nu_1^{(2)},\dots,\nu_n^{(2)}]],\dots, [\nu_1^{(S)},\dots,\nu_n^{(S)}]]\)
    type(mol_block)         :: res
    type(group_permutation) :: s
    integer(IK)             :: napm
!
    s = group_permutation(sym)
    napm = MAX(n, 0)
    if (PRESENT(sym)) napm = MAX(napm, SIZE(sym, 1))
    res%q = [napm, MAX(M, 1), s%q]
!
  end function mol_block_new
!
!| number of atoms, defined by \(nM\).
  pure function mol_block_natm(q) result(res)
    integer(IK), intent(in)     :: q(*)
    !! mol_block
    integer(IK)                 :: res
    res = q(pm) * q(pn)
  end function mol_block_natm
!
!| number of molecules, \(M\).
  pure function mol_block_nmol(q) result(res)
    integer(IK), intent(in)     :: q(*)
    !! mol_block
    integer(IK)                 :: res
    res = q(pm)
  end function mol_block_nmol
!
!| number of atoms per molecule, \(n\)
  pure function mol_block_napm(q) result(res)
    integer(IK), intent(in)     :: q(*)
    !! mol_block
    integer(IK)                 :: res
    res = q(pn)
  end function mol_block_napm
!
!| number of molecular symmetry, \(S\).
  pure function mol_block_nsym(q) result(res)
    integer(IK), intent(in)     :: q(*)
    !! mol_block
    integer(IK)                 :: res
    res = group_permutation_nsym(q(pg))
  end function mol_block_nsym
!
!| memory blocksize per molecule, defined by \(dn\).
  pure function mol_block_each_size(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! mol_block
    integer(IK)             :: res
    res = D * q(pn)
  end function mol_block_each_size
!
!| memory blocksize, defined by \(dnM\).
  pure function mol_block_total_size(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! mol_block
    integer(IK)             :: res
    res = D * q(pn) * q(pm)
  end function mol_block_total_size
!
!| Compute molecular symmetry permutation according to isym.
  pure subroutine mol_block_swap(q, isym, X)
    integer(IK), intent(in) :: q(*)
    !! mol_block
    integer(IK), intent(in) :: isym
    !! symmetry id
    real(RK), intent(inout) :: X(*)
    !! work array
!
    call group_permutation_swap(q(pg), isym, D, X)
!
  end subroutine mol_block_swap
!
!| Compute molecular symmetry inverse permutation according to isym.
  pure subroutine mol_block_inverse_swap(q, isym, X)
    integer(IK), intent(in) :: q(*)
    !! mol_block
    integer(IK), intent(in) :: isym
    !! symmetry id.
    real(RK), intent(inout) :: X(*)
    !! work array
!
    call group_permutation_inverse(q(pg), isym, D, X)
!
  end subroutine mol_block_inverse_swap
!
end module mod_mol_block

