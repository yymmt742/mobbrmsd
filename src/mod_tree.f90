module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_molecular_rotation
  use mod_mol_block
  use mod_partial_rmsd
  implicit none
  private
  public :: node, breadth
!
  type node
    private
    type(partial_rmsd)             :: prd
    type(mol_block_list)           :: blk
    integer(IK), allocatable       :: per(:)
    !! molecular permutation block
    real(RK), public               :: lowerbound = RHUGE
    !! the lower bound
  contains
!   procedure         :: generate_breadth => node_generate_breadth
    final             :: node_destroy
  end type node
!
  type breadth
    type(node), allocatable :: nodes(:)
  contains
    final             :: breadth_destroy
  end type breadth
!
  interface node
    module procedure node_new
  end interface node
!
contains
!
! pure elemental function nwork(blk) result(res)
!   class(mol_block_list), intent(in) :: blk
!   integer(IK)                       :: ispc, res
!   ispc = blk%ispecies()
!   res = 2 * blk%d * blk%b(ispc)%m
! end function nwork
!
!| generate node instance
  pure function node_new(blk, x, y) result(res)
    class(mol_block_list), intent(in) :: blk
    real(RK), intent(in)              :: X(*)
    real(RK), intent(in)              :: Y(*)
    type(node)                        :: res
!
    res%blk = blk
    call proc_fixed_part(blk, X, Y, res%prd)
    res%lowerbound = res%prd%sd()
!
  end function node_new
!
  pure subroutine proc_fixed_part(blk, X, Y, prd)
    class(mol_block_list), intent(in) :: blk
    real(RK), intent(in)              :: X(*)
    real(RK), intent(in)              :: Y(*)
    type(partial_rmsd), intent(inout) :: prd
    real(RK)                          :: W(nwork1(blk))
    integer(IK)                       :: r(blk%nspecies())
    integer(IK)                       :: p(blk%nspecies())
    integer(IK)                       :: q(blk%nspecies())
    integer(IK)                       :: i, ix, iy, s, n
!
    s = blk%nspecies()
!
    do i = 1, s
      p(i) = blk%res_pointer(i)
    end do
!
    r = blk%n_res()
    n = SUM(r)
    r = blk%d * r
!
    q(1) = 0
    do i = 1, s - 1
      q(i + 1) = r(i) + q(i)
    end do
!
    ix = 1
    iy = ix + blk%d * n
!
    do concurrent(i = 1:s)
      if (r(i) > 0) call copy(r(i), X(p(i)), w(ix + q(i)))
    end do
    do concurrent(i = 1:s)
      if (r(i) > 0) call copy(r(i), Y(p(i)), w(iy + q(i)))
    end do
    prd = partial_rmsd(blk%d, n, W(ix:), W(iy:))
  end subroutine proc_fixed_part
!
  pure elemental function nwork1(blk) result(res)
    class(mol_block_list), intent(in) :: blk
    integer(IK)                       :: res
    res = MAX(1, 2 * SUM(blk%n_res()) * blk%d)
  end function nwork1
!
  pure subroutine copy(n, x, w)
    integer(IK), intent(in) :: n
    real(RK), intent(in)    :: x(*)
    real(RK), intent(inout) :: w(*)
    integer(IK)             :: i
    do concurrent(i=1:n)
      w(i) = x(i)
    end do
  end subroutine copy
!
! pure subroutine centering(d, n, x)
!   integer(IK), intent(in) :: d, n
!   real(RK), intent(inout) :: x(d, n)
!   real(RK)                :: c(d)
!   integer(IK)             :: i, j
!   c = SUM(X, 2) / n
!   do concurrent(i=1:d, j=1:n)
!     x(i, j) = x(i, j) - c(i)
!   end do
! end subroutine centering
!
!| generate childe nodes instance
!  function node_generate_breadth(this, rot, x, y) result(res)
!   class(node), intent(in)               :: this
!   !! molecular template
!   class(molecular_rotation), intent(in) :: rot(:)
!   !! molecular template
!   real(RK), intent(in)                  :: x(*)
!   !! reference molecular coordinate, x(d,m,n)
!   real(RK), intent(in)                  :: y(*)
!   !! target molecular coordinate, y(d,m,n)
!   type(breadth)                         :: res
!   type(mol_block_list)                  :: b_child
!   integer(IK)                           :: nper, nsym, nspc
!   integer(IK)                           :: iper, isym, ispc
!
!   b_child = this%blk%child()
!   ispc = b_child%ispecies()
!   nspc = b_child%nspecies()
!   if (ispc < 1) ispc = nspc
!   nper = b_child%b(ispc)%g + 1
!   nsym = rot(ispc)%n_sym()
!
!   allocate (res%nodes(nper * (nsym + 1)))
!
!   do isym=0,nsym
!   do iper=1,nper
!   !do concurrent(iper=1:nper, isym=0:nsym)
!     block
!       type(molecular_permutation) :: p_child(nspc)
!       integer(IK)                 :: i
!       do concurrent(i=1:ispc - 1)
!         p_child(i) = this%per(i)
!       end do
!       p_child(ispc) = this%per(ispc)%child(iper, isym)
!       do concurrent(i=ispc + 1:nspc)
!         p_child(i) = this%per(i)
!       end do
!       res%nodes(isym * nper + iper) = node(b_child, rot, p_child, x, y)
!     end block
!   end do
!   end do
!
! end function node_generate_breadth
!
  pure elemental subroutine node_destroy(this)
    type(node), intent(inout) :: this
    call this%prd%clear()
    if (ALLOCATED(this%per)) deallocate (this%per)
    this%lowerbound = ZERO
  end subroutine node_destroy
!
  pure elemental subroutine breadth_destroy(this)
    type(breadth), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
  end subroutine breadth_destroy
!
end module mod_tree
