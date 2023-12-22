module mod_branch_and_prune
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_mol_block
  use mod_molecular_rotation
  use mod_d_matrix
  use mod_tree
  implicit none
  private
  public :: branch_and_prune
!
  type branch_and_prune
    type(d_matrix_list)        :: dmat
    type(molecular_rotation), allocatable :: rot(:)
    type(node)                 :: root
    type(breadth), allocatable :: childs(:)
    real(RK)                   :: upperbound = RHUGE
  contains
    procedure :: memsize    => branch_and_prune_memsize
    procedure :: setup      => branch_and_prune_setup
    procedure :: run        => branch_and_prune_run
    procedure :: clear      => branch_and_prune_clear
    final     :: branch_and_prune_destroy
  end type branch_and_prune
!
  interface branch_and_prune
    module procedure branch_and_prune_new
  end interface branch_and_prune
!
contains
!
!| generate node instance
  pure function branch_and_prune_new(blk, rot, p) result(res)
    type(mol_block_list), intent(in)     :: blk
    type(molecular_rotation), intent(in) :: rot(*)
    integer(IK), intent(in)              :: p
    type(branch_and_prune)               :: res
    integer(IK)                          :: n, i
!
    res%dmat = d_matrix_list(blk, p)
    n = res%dmat%n_depth()
    allocate (res%childs(n))
!
    n = SIZE(blk%b)
    allocate (res%rot(n))
!
    do concurrent(i=1:n)
      res%rot(i) = rot(i)
    end do
!
  end function branch_and_prune_new
!
  pure elemental function branch_and_prune_memsize(this) result(res)
    class(branch_and_prune), intent(in) :: this
    integer(IK)                         :: res
    res = this%dmat%memsize()
  end function branch_and_prune_memsize
!
  pure subroutine branch_and_prune_setup(this, X, Y, W)
    class(branch_and_prune), intent(inout) :: this
    real(RK), intent(in)                   :: X(*)
    real(RK), intent(in)                   :: Y(*)
    real(RK), intent(inout)                :: W(*)
!
    call this%dmat%eval(this%rot, X, Y, W)
    this%root = node(this%dmat, W)
    this%childs(1) = this%root%generate_breadth(this%dmat, W)
    this%upperbound = RHUGE
!
  end subroutine branch_and_prune_setup
!
   subroutine branch_and_prune_run(this, W)
    class(branch_and_prune), intent(inout) :: this
    real(RK), intent(in)                   :: W(*)
    integer(IK)                            :: ndep
!
    if (.not. ALLOCATED(this%childs)) return
    if (SIZE(this%childs) < 1) return
!
    ndep = SIZE(this%childs)
    block
      integer(IK) :: bstper(ndep)
      call breadth_search(1, ndep, this%dmat, W, this%childs, this%upperbound, bstper)
    end block
!
  end subroutine branch_and_prune_run
!
  recursive subroutine breadth_search(idep, ndep, dmat, W, childs, upperbound, bstper)
    integer(IK), intent(in)         :: idep, ndep
    type(d_matrix_list), intent(in) :: dmat
    real(RK), intent(in)            :: W(*)
    type(breadth), intent(inout)    :: childs(*)
    real(RK), intent(inout)         :: upperbound
    integer(IK), intent(inout)      :: bstper(*)
    integer(IK)                     :: i
!
print*,childs(idep)%is_finished(),idep, childs(idep)%nodes%alive
print'(*(f9.3))', upperbound, childs(idep)%nodes%lowerbound
    if (idep == ndep) then
      block
        real(RK) :: lb
        call childs(idep)%set_node_minloc()
        lb = childs(idep)%lowerbound()
        if (upperbound > lb) then
          upperbound = lb
          do i = 1, ndep
            print *, childs(i)%iper(), childs(i)%isym()
          end do
        end if
      end block
      return
    end if
!
    do while (childs(idep)%not_finished())
      call childs(idep)%set_node_index()
      childs(idep + 1) = childs(idep)%generate_breadth(dmat, W)
      call breadth_search(idep + 1, ndep, dmat, W, childs, upperbound, bstper)
      call childs(idep)%prune(upperbound)
      print*,childs(idep)%nodes%alive
    end do
!
  end subroutine breadth_search
!
  pure elemental subroutine branch_and_prune_clear(this)
    class(branch_and_prune), intent(inout) :: this
    call this%dmat%clear()
    call this%root%clear()
    if (ALLOCATED(this%childs)) deallocate (this%childs)
    if (ALLOCATED(this%rot)) deallocate (this%rot)
    this%upperbound = RHUGE
  end subroutine branch_and_prune_clear
!
  pure elemental subroutine branch_and_prune_destroy(this)
    type(branch_and_prune), intent(inout) :: this
    call branch_and_prune_clear(this)
  end subroutine branch_and_prune_destroy
!
end module mod_branch_and_prune
