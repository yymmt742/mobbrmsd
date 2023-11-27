module mod_molecular_topology
  use mod_params, only: IK, RK
  use mod_base
  use mod_element
  use mod_element_set
  use mod_molecule
  implicit none
  private
  public :: molecular_topology
!
  type, extends(base) :: molecular_topology
    private
    type(element_set), allocatable :: e
    type(molecule), allocatable    :: m(:)
    integer(IK), allocatable       :: l(:)
  contains
    procedure :: nmol       => molecular_topology_nmol
    procedure :: nspecies   => molecular_topology_nspecies
    procedure :: clear      => molecular_topology_clear
    final     :: molecular_topology_destroy
  end type molecular_topology
!
  interface molecular_topology
    module procedure molecular_topology_new
  end interface molecular_topology
!
contains
!
  pure function molecular_topology_new(ele, mol, list) result(res)
    class(element_set), intent(in) :: ele
    class(molecule), intent(in)    :: mol(:)
    integer(IK), intent(in)        :: list(:)
    type(molecular_topology)       :: res
    integer(IK)                    :: i, s, n
!
    s = SIZE(mol)
    n = SIZE(list)
    allocate (res%m(s))
    do concurrent(i=1:s)
      res%m(i) = mol(i)
    end do
!
    allocate (res%l, source=list)
!
  end function molecular_topology_new
!
  pure elemental function molecular_topology_nmol(this) result(res)
    class(molecular_topology), intent(in) :: this
    integer(IK)                     :: res
    if (ALLOCATED(this%l)) then
      res = SIZE(this%l)
    else
      res = 0
    end if
  end function molecular_topology_nmol
!
  pure elemental function molecular_topology_nspecies(this) result(res)
    class(molecular_topology), intent(in) :: this
    integer(IK)                     :: res
    if (ALLOCATED(this%m)) then
      res = SIZE(this%m)
    else
      res = 0
    end if
  end function molecular_topology_nspecies
!
  pure elemental subroutine molecular_topology_clear(this)
    class(molecular_topology), intent(inout) :: this
    call this%e%clear()
    if (ALLOCATED(this%m)) deallocate (this%m)
    if (ALLOCATED(this%l)) deallocate (this%l)
  end subroutine molecular_topology_clear
!
  pure elemental subroutine molecular_topology_destroy(this)
    type(molecular_topology), intent(inout) :: this
    call this%clear()
  end subroutine molecular_topology_destroy
!
end module mod_molecular_topology
