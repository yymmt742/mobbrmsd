module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_molecular_rotation
  use mod_mol_block
  use mod_Hungarian
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
    procedure         :: generate_breadth => node_generate_breadth
    procedure         :: has_child        => node_has_child
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
  pure subroutine calc_lowerbound_matrix(blk, X, Y, Z)
    class(mol_block_list), intent(in)     :: blk
    real(RK), intent(in)                  :: X(*)
    real(RK), intent(in)                  :: Y(*)
    real(RK), intent(inout)               :: Z(*)
    integer(IK)                           :: si
    integer(IK)                           :: nspc, ispc
!
    nspc = blk%nspecies()
    do concurrent(ispc=1:nspc)
      block
        integer(IK) :: s, px, pz
        s = blk%r(ispc)%n_sym() + 1
        px = bx(ispc, blk)
        pz = bz(ispc, blk)
        ph = pz + blk%b(ispc)%g**2 * s
        call lbmat(blk%d, s, blk%b(ispc), blk%r(ispc), X(px), Y(px), Z(pz))
      end block
    end do
!
  contains
!
    pure elemental function bx(s, blk) result(res)
      integer(IK), intent(in)           :: s
      class(mol_block_list), intent(in) :: blk
      integer(IK)                       :: i, res
      res = 0
      do i = 1, s - 1
        res = res + blk%b(i)%m * blk%b(i)%n
      end do
      res = res * blk%d + 1
    end function bx
!
    pure elemental function bz(s, blk) result(res)
      integer(IK), intent(in)           :: s
      class(mol_block_list), intent(in) :: blk
      integer(IK)                       :: i, res
      res = s + s - 1
      do i = 1, s - 1
        res = res + blk%b(i)%g**2
      enddo
      do i = 1, s - 1
        res = res + blk%r(i)%n_sym()
      enddo
    end function bz
!
  end subroutine calc_lowerbound_matrix
!
  pure subroutine lbmat(d, s, b, r, X, Y, Z)
    integer(IK), intent(in)              :: d, s
    type(mol_block), intent(in)          :: b
    type(molecular_rotation), intent(in) :: r
    real(RK), intent(in)                 :: X(d, b%m, b%g), Y(d, b%m, b%g)
    real(RK), intent(inout)              :: Z(b%g, b%g, s)
    integer(IK)                          :: i, j, k, dm
!
    dm = d * m
    do concurrent(i=1:b%g, j=1:b%g, k=1:s)
      if (i < j) CYCLE
      block
        type(partial_rmsd) :: prd
        real(RK)           :: W(dm, 2)
        call copy(dm, X(1, 1, i), W(1, 1))
        call copy(dm, Y(1, 1, i), W(1, 2))
        if (k > 1) call r%swap(d, W(1, 2), k - 1)
        prd = partial_rmsd(d, m, W(1, 1), W(1, 2))
        Z(i, j, k) = prd%sd()
        if (i > j) Z(j, i, k) = Z(i, j, k)
      end block
    enddo
!
  end subroutine lbmat
!
!| generate node instance
  pure function node_new(blk, x, y) result(res)
    class(mol_block_list), intent(in) :: blk
    real(RK), intent(in)              :: X(*)
    real(RK), intent(in)              :: Y(*)
    type(node)                        :: res
    integer(IK)                       :: si, nb
!
    res%blk = blk
    si = blk%ispecies()
    nb = 0; if (si > 0) nb = blk%b(si)%g
    allocate (res%per(0))
!
    call proc_fixed_part(blk, X, Y, res%prd)
    res%lowerbound = res%prd%sd()
!
  end function node_new
!
  pure elemental function node_has_child(this) result(res)
    class(node), intent(in) :: this
    logical                 :: res
    res = this%blk%has_child()
  end function node_has_child
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
!
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
!| generate childe nodes instance
  function node_generate_breadth(this, rot, x, y) result(res)
    class(node), intent(in)               :: this
    !! molecular template
    class(molecular_rotation), intent(in) :: rot(:)
    !! molecular template
    real(RK), intent(in)                  :: x(*)
    !! reference molecular coordinate, x(d,m,n)
    real(RK), intent(in)                  :: y(*)
    !! target molecular coordinate, y(d,m,n)
    type(breadth)                         :: res
    type(mol_block_list)                  :: b_child
    integer(IK), allocatable              :: per(:)
    integer(IK)                           :: px, bm, bs
    integer(IK)                           :: iper, ipnt, isym, ispc
    integer(IK)                           :: nper, nsym
!
    if (.not. this%has_child()) then
      allocate (res%nodes(0))
      return
    end if
!
    b_child = this%blk%child()
    ispc = this%blk%ispecies()
!
    ipnt = this%blk%ipointer(0)
    nper = this%blk%b(ispc)%g
    nsym = rot(ispc)%n_sym()
!
    bm = this%blk%b(ispc)%m
    bs = this%blk%d * bm
    px = ipnt + SIZE(this%per) * bs
!
    per = resper(nper, nper + SIZE(this%per), this%per)
    allocate (res%nodes(nper * (nsym + 1)))
!
!print*,ipnt, nper, nsym, bm, bs, px
print*, per
!
    do concurrent(iper=1:nper, isym=0:nsym)
      block
        integer(IK) :: py, inod
        real(RK)    :: w(bs, 2)
        inod = isym * nper + iper
        res%nodes(inod)%blk = b_child
        if (nper == 1) then
          allocate (res%nodes(inod)%per(0))
        else
          res%nodes(inod)%per = [this%per, per(iper)]
        end if
        py = ipnt + (per(iper) - 1) * bs
        call copy(bs, X(px), w(1, 1))
        call copy(bs, Y(py), w(1, 2))
        call rot(ispc)%swap(b_child%d, w(1, 2), isym)
        res%nodes(inod)%prd = this%prd%append(bm, w(:, 1), w(:, 2))
        res%nodes(inod)%lowerbound = res%nodes(inod)%prd%sd()
      end block
    end do
!
  end function node_generate_breadth
!
  pure function resper(g, n, per) result(res)
    integer(IK), intent(in) :: g, n, per(n)
    integer(IK)             :: i, j, res(g)
    j = 0
    do i = 1, n
      if (ANY(i == per)) cycle
      j = j + 1
      res(j) = i
      if (j == g) return
    end do
  end function resper
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
