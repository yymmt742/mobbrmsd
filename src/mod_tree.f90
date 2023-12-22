module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_molecular_rotation
  use mod_mol_block
  use mod_d_matrix
  implicit none
  private
  public :: node, breadth
!
  type node
    private
    integer(IK)              :: depth = 0
    integer(IK)              :: ispc = 0
    integer(IK), allocatable :: per(:)
    real(RK)                 :: H
    real(RK), allocatable    :: C(:), L(:)
    real(RK), public         :: lowerbound = ZERO
    logical, public          :: alive = .TRUE.
    !! the lower bound
  contains
    procedure         :: n_per            => node_n_per
    procedure         :: n_sym            => node_n_sym
    procedure         :: n_depth          => node_n_depth
    procedure         :: generate_breadth => node_generate_breadth
    procedure         :: has_child        => node_has_child
    procedure         :: clear            => node_clear
    final             :: node_destroy
  end type node
!
  type breadth
    type(node), allocatable :: nodes(:)
  contains
    procedure         :: prune            => breadth_prune
    procedure         :: is_finished      => breadth_is_finished
    procedure         :: lowervalue_index => breadth_lowervalue_index
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
!
    res%L = res%L(res%ispc + 1:)
    res%alive = res%has_child()
!
  end function node_new_as_root
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
    res%alive = res%has_child()
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
  pure elemental function node_n_depth(this, dmat) result(res)
    class(node), intent(in)         :: this
    type(d_matrix_list), intent(in) :: dmat
    integer(IK)                     :: res
    if(this%ispc<1.or.dmat%l<this%ispc)then
      res = 0
    else
      res = SUM(dmat%m(this%ispc:)%g) - this%depth + 1
    endif
  end function node_n_depth
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
  pure elemental function node_has_child(this) result(res)
    class(node), intent(in) :: this
    logical                 :: res
    res = this%ispc > 0
  end function node_has_child
!
!| generate childe nodes instance
  pure function node_generate_breadth(this, dmat, W) result(res)
    class(node), intent(in)         :: this
    type(d_matrix_list), intent(in) :: dmat
    real(RK), intent(in)            :: W(*)
    type(breadth)                   :: res
    integer(IK)                     :: nper, nsym, nchd, iper, isym
!
    if (.not. this%has_child()) then
      allocate (res%nodes(0))
      return
    end if
!
    nper = this%n_per(dmat)
    nsym = this%n_sym(dmat)
    nchd = nper * nsym
!
    allocate (res%nodes(nper * nsym))
!
    do concurrent(iper=1:nper, isym=1:nsym)
      block
        integer(IK) :: inod
        inod = (isym - 1) * nper + iper
        res%nodes(inod) = node_new(this, iper, isym, dmat, W)
      end block
    end do
!
  end function node_generate_breadth
!
  pure elemental subroutine node_clear(this)
    class(node), intent(inout) :: this
    if (ALLOCATED(this%per)) deallocate (this%per)
    if (ALLOCATED(this%L)) deallocate (this%L)
    if (ALLOCATED(this%C)) deallocate (this%C)
    this%depth = 0
    this%ispc  = 0
    this%lowerbound = ZERO
  end subroutine node_clear
!
  pure elemental subroutine node_destroy(this)
    type(node), intent(inout) :: this
    call node_clear(this)
  end subroutine node_destroy
!
  pure elemental subroutine breadth_prune(this, upperbound)
    class(breadth), intent(inout) :: this
    real(RK), intent(in)          :: upperbound
    integer(IK)                   :: i, n
    if (ALLOCATED(this%nodes)) then
      n = SIZE(this%nodes)
      do concurrent(i=1:n)
        this%nodes(i)%alive = this%nodes(i)%alive .and. this%nodes(i)%lowerbound < upperbound
      end do
    end if
  end subroutine breadth_prune
!
  pure elemental function breadth_is_finished(this) result(res)
    class(breadth), intent(in) :: this
    logical                    :: res
    if (ALLOCATED(this%nodes)) then
      res = .not. ANY(this%nodes%alive)
    else
      res = .true.
    end if
  end function breadth_is_finished
!
  pure elemental function breadth_lowervalue_index(this) result(res)
    class(breadth), intent(in) :: this
    real(RK)                   :: lowervalue
    integer(IK)                :: res, i, n
!
    res = 0
    if (.not.ALLOCATED(this%nodes)) return
!
    lowervalue = RHUGE
    n = SIZE(this%nodes)
    do i = 1, n
      if (this%nodes(i)%alive .and. this%nodes(i)%lowerbound < lowervalue) then
        res = i
        lowervalue = this%nodes(i)%lowerbound
      end if
    end do
!
  end function breadth_lowervalue_index

  pure elemental subroutine breadth_clear(this)
    class(breadth), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
  end subroutine breadth_clear
!
  pure elemental subroutine breadth_destroy(this)
    type(breadth), intent(inout) :: this
    call breadth_clear(this)
  end subroutine breadth_destroy
!
end module mod_tree
