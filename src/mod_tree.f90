module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_molecular_rotation
  use mod_mol_block
  use mod_d_matrix
  use mod_partial_rmsd
  implicit none
  private
  public :: d_matrix_list, node, breadth
!
  type node
    private
    integer(IK)              :: depth = 0
    integer(IK)              :: ispc = 0
    integer(IK), allocatable :: per(:)
    real(RK)                 :: H
    real(RK), allocatable    :: C(:), L(:)
    !! molecular permutation block
    real(RK), public         :: lowerbound = RHUGE
    !! the lower bound
  contains
    procedure         :: n_per            => node_n_per
    procedure         :: n_sym            => node_n_sym
!   procedure         :: generate_breadth => node_generate_breadth
!   procedure         :: has_child        => node_has_child
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
    module procedure node_new, node_new_as_root
  end interface node
!
contains
!| generate node instance
  pure function node_new_as_root(dmat, W) result(res)
    type(d_matrix_list), intent(in) :: dmat
    real(RK), intent(in)            :: W(*)
    type(node)                      :: res
    integer(IK)                     :: i, g
!
    res%H = W(dmat%h)
!
    allocate (res%C, source=W(dmat%c:dmat%c + dmat%d**2 - 1))
    allocate (res%L(dmat%l))
!
    do concurrent(i=1:dmat%l)
      block
        real(RK)    :: LF, LB, H, C(1)
        integer(IK) :: ires(dmat%m(i)%g)
        integer(IK) :: j
        do concurrent(j=1:dmat%m(i)%g)
          ires(j) = j
        end do
        call d_matrix_partial_eval(dmat%m(i), 0, 0, 0, ires, W, LF, LB, H, C)
        res%L(i) = LB
      end block
    end do
!
    res%lowerbound = W(dmat%v) + SUM(res%L)
!
    call next_index(dmat%l, dmat%m, 0, 0, res)
    if(res%ispc<1) return
!
    g = dmat%m(res%ispc)%g
    allocate (res%per(g))
    do concurrent(i=1:g)
      res%per(i) = i
    end do
    res%L = res%L(res%ispc + 1:)
!
  end function node_new_as_root
!
  pure subroutine next_index(l, m, depth, ispc, child)
    integer(IK), intent(in)    :: l, depth, ispc
    type(d_matrix), intent(in) :: m(l)
    type(node), intent(inout)  :: child
    integer(IK)                :: i, spc
    spc = MAX(ispc, 1)
    if (depth < m(spc)%g) then
      child%ispc = spc
      child%depth = depth + 1
      return
    end if
    do i = spc + 1, l
      if (m(i)%g < 1) cycle
      child%ispc = i
      child%depth = 1
      return
    end do
  end subroutine next_index
!
  pure elemental function node_n_sym(this, dmat) result(res)
    class(node), intent(in)         :: this
    type(d_matrix_list), intent(in) :: dmat
    integer(IK)                     :: res
    if(this%ispc<1.or.dmat%l<this%ispc)then
      res = 0
    else
      res = dmat%m(this%ispc)%s
    endif
  end function node_n_sym
!
  pure elemental function node_n_per(this, dmat) result(res)
    class(node), intent(in)         :: this
    type(d_matrix_list), intent(in) :: dmat
    integer(IK)                     :: res
    if (this%ispc < 1 .or. dmat%l < this%ispc) then
      res = 0
    else
      res = dmat%m(this%ispc)%g - this%depth + 1
    end if
  end function node_n_per
!
  pure function node_new(parent, iper, isym, dmat, W) result(res)
    type(node), intent(in)          :: parent
    type(d_matrix_list), intent(in) :: dmat
    integer(IK), intent(in)         :: iper, isym
    real(RK), intent(in)            :: W(*)
    real(RK)                        :: LF, LB
    type(node)                      :: res
    integer(IK)                     :: i, np
!
    allocate (res%C, source=parent%C)
    res%H = parent%H
!
    np = SIZE(parent%per) - parent%depth
    block
      integer(IK) :: ip, ir(np)
!
      call perm_index(np, iper, parent%per(parent%depth), ip, ir)
      call d_matrix_partial_eval(dmat%m(parent%ispc), parent%depth, ip, isym, ir, &
     &                           W, LF, LB, res%H, res%C)
!
      res%lowerbound = LF + LB + SUM(parent%L)
!
      call next_index(dmat%l, dmat%m, parent%depth, parent%ispc, res)
!
      if (parent%ispc == res%ispc) then
        ALLOCATE(res%per, mold=parent%per)
        do i = 1, parent%depth - 1
          res%per(i) = parent%per(i)
        end do
        res%per(parent%depth) = ip
        do i = parent%depth + 1, SIZE(res%per)
          res%per(i) = ir(i - parent%depth)
        enddo
        allocate (res%L, source=parent%L)
      elseif (parent%ispc < res%ispc) then
        allocate (res%per(dmat%m(res%ispc)%g))
        do concurrent(i=1:dmat%m(res%ispc)%g)
          res%per(i) = i
        end do
        i = res%ispc - parent%ispc + 1
        allocate (res%L, source=parent%L(i:))
      else
        allocate (res%L(0))
      end if
    end block
!
  contains
!
    pure subroutine perm_index(n, j, q, p, r)
      integer(IK), intent(in)    :: n, j, q(n + 1)
      integer(IK), intent(inout) :: p, r(n)
      integer(IK)                :: i
      p = q(j)
      do concurrent(i=1:j-1)
        r(i) = q(i)
      end do
      do concurrent(i=j:n)
        r(i) = q(i + 1)
      end do
    end subroutine perm_index
!
  end function node_new
!
! pure elemental function node_has_child(this) result(res)
!   class(node), intent(in) :: this
!   logical                 :: res
!   res = this%blk%has_child()
! end function node_has_child
!
! pure subroutine proc_fixed_part(blk, X, Y, prd)
!   class(mol_block_list), intent(in) :: blk
!   real(RK), intent(in)              :: X(*)
!   real(RK), intent(in)              :: Y(*)
!   type(partial_rmsd), intent(inout) :: prd
!   real(RK)                          :: W(nwork1(blk))
!   integer(IK)                       :: r(blk%nspecies())
!   integer(IK)                       :: p(blk%nspecies())
!   integer(IK)                       :: q(blk%nspecies())
!   integer(IK)                       :: i, ix, iy, s, n
!
!   s = blk%nspecies()
!
!   do i = 1, s
!     p(i) = blk%res_pointer(i)
!   end do
!
!   r = blk%n_res()
!   n = SUM(r)
!   r = blk%d * r
!
!   q(1) = 0
!   do i = 1, s - 1
!     q(i + 1) = r(i) + q(i)
!   end do
!
!   ix = 1
!   iy = ix + blk%d * n
!
!   do concurrent(i = 1:s)
!     if (r(i) > 0) call copy(r(i), X(p(i)), w(ix + q(i)))
!   end do
!   do concurrent(i = 1:s)
!     if (r(i) > 0) call copy(r(i), Y(p(i)), w(iy + q(i)))
!   end do
!   prd = partial_rmsd(blk%d, n, W(ix:), W(iy:))
!
! end subroutine proc_fixed_part
!
! pure elemental function nwork1(blk) result(res)
!   class(mol_block_list), intent(in) :: blk
!   integer(IK)                       :: res
!   res = MAX(1, 2 * SUM(blk%n_res()) * blk%d)
! end function nwork1
!
! pure subroutine copy(n, x, w)
!   integer(IK), intent(in) :: n
!   real(RK), intent(in)    :: x(*)
!   real(RK), intent(inout) :: w(*)
!   integer(IK)             :: i
!   do concurrent(i=1:n)
!     w(i) = x(i)
!   end do
! end subroutine copy
!
!| generate childe nodes instance
! function node_generate_breadth(this, rot, x, y) result(res)
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
!   integer(IK), allocatable              :: per(:)
!   integer(IK)                           :: px, bm, bs
!   integer(IK)                           :: iper, ipnt, isym, ispc
!   integer(IK)                           :: nper, nsym
!
!   if (.not. this%has_child()) then
!     allocate (res%nodes(0))
!     return
!   end if
!
!   b_child = this%blk%child()
!   ispc = this%blk%ispecies()
!
!   ipnt = this%blk%ipointer(0)
!   nper = this%blk%b(ispc)%g
!   nsym = rot(ispc)%n_sym()
!
!   bm = this%blk%b(ispc)%m
!   bs = this%blk%d * bm
!   px = ipnt + SIZE(this%per) * bs
!
!   per = resper(nper, nper + SIZE(this%per), this%per)
!   allocate (res%nodes(nper * (nsym + 1)))
!
!print*,ipnt, nper, nsym, bm, bs, px
!rint*, per
!
!   do concurrent(iper=1:nper, isym=0:nsym)
!     block
!       integer(IK) :: py, inod
!       real(RK)    :: w(bs, 2)
!       inod = isym * nper + iper
!       res%nodes(inod)%blk = b_child
!       if (nper == 1) then
!         allocate (res%nodes(inod)%per(0))
!       else
!         res%nodes(inod)%per = [this%per, per(iper)]
!       end if
!       py = ipnt + (per(iper) - 1) * bs
!       call copy(bs, X(px), w(1, 1))
!       call copy(bs, Y(py), w(1, 2))
!       call rot(ispc)%swap(b_child%d, w(1, 2), isym)
!       res%nodes(inod)%prd = this%prd%append(bm, w(:, 1), w(:, 2))
!       res%nodes(inod)%lowerbound = res%nodes(inod)%prd%sd()
!     end block
!   end do
!
! end function node_generate_breadth
!
! pure function resper(g, n, per) result(res)
!   integer(IK), intent(in) :: g, n, per(n)
!   integer(IK)             :: i, j, res(g)
!   j = 0
!   do i = 1, n
!     if (ANY(i == per)) cycle
!     j = j + 1
!     res(j) = i
!     if (j == g) return
!   end do
! end function resper
!
  pure elemental subroutine node_destroy(this)
    type(node), intent(inout) :: this
    if (ALLOCATED(this%per)) deallocate (this%per)
    if (ALLOCATED(this%L)) deallocate (this%L)
    if (ALLOCATED(this%C)) deallocate (this%C)
    this%depth = 0
    this%ispc  = 0
    this%lowerbound = ZERO
  end subroutine node_destroy
!
  pure elemental subroutine breadth_destroy(this)
    type(breadth), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
  end subroutine breadth_destroy
!
end module mod_tree
