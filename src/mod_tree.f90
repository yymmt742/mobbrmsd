module mod_tree
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_molecular_rotation
  use mod_mol_block
  use mod_d_matrix
  implicit none
  private
  public :: tree
  public :: node, breadth
!
  type node_
    private
    sequence
    logical     :: alive
    integer(IK) :: iper, isym, v, h, c
  end type node_
!
  type breadth_
    private
    integer(IK)              :: inod = 0
    integer(IK)              :: nper = 0
    integer(IK)              :: nsym = 0
    integer(IK)              :: nnod = 0
    integer(IK)              :: ispc = 0
    integer(IK)              :: irow = 0
    integer(IK)              :: lowd = 0
    integer(IK)              :: uppd = 0
    type(node_), allocatable :: nodes(:)
  contains
    final             :: breadth__destroy
  end type breadth_
!
  type tree
    private
    integer(IK)                 :: iscope
    type(d_matrix_list)         :: dmat
    type(node_)                 :: root
    integer(IK), allocatable    :: perms(:)
    type(breadth_), allocatable :: breadthes(:)
    type(molecular_rotation), allocatable :: rots(:)
  contains
    procedure         :: init             => tree_init
    procedure         :: setup            => tree_setup
    procedure         :: n_depth          => tree_n_depth
    procedure         :: memsize          => tree_memsize
    final             :: tree_destroy
  end type tree
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
    integer(IK)             :: inod = 0
    integer(IK)             :: nper = 0
    integer(IK)             :: nsym = 0
    integer(IK)             :: nnod = 0
    type(node), allocatable :: nodes(:)
  contains
    procedure         :: generate_breadth => breadth_generate_breadth
    procedure         :: prune            => breadth_prune
    procedure         :: is_finished      => breadth_is_finished
    procedure         :: not_finished     => breadth_not_finished
    procedure         :: isym             => breadth_isym
    procedure         :: iper             => breadth_iper
    procedure         :: set_node_index   => breadth_set_node_index
    procedure         :: set_node_minloc  => breadth_set_node_minloc
    procedure         :: lowerbound       => breadth_lowerbound
    final             :: breadth_destroy
  end type breadth
!
  interface node
    module procedure node_new, node_new_as_root
  end interface node
!
contains
!
  pure subroutine node__init_as_root(this, dmat, pw)
    type(node_), intent(inout)      :: this
    type(d_matrix_list), intent(in) :: dmat
    integer(IK), intent(in)         :: pw
!
    this%iper = 0
    this%isym = 0
    this%v = pw          ! v(1)
    this%h = this%v + 1  ! h(1)
    this%c = this%h + 1  ! c(d, d)
!
  end subroutine node__init_as_root
!
  pure subroutine node__init(this, iper, isym, pw)
    type(node_), intent(inout)      :: this
    integer(IK), intent(in)         :: iper, isym, pw
!
    this%iper = iper
    this%isym = isym
    this%v = pw          ! v(1)
    this%h = this%v + 1  ! h(1)
    this%c = this%h + 1  ! c(d, d)
!
  end subroutine node__init
!
  pure elemental function node__memsize(dmat) result(res)
    type(d_matrix_list), intent(in) :: dmat
    integer(IK)                     :: res
    res = 2 + dmat%d**2
  end function node__memsize
!
  pure subroutine breadth__init(this, dmat, pw, depth)
    type(breadth_), intent(inout)   :: this
    type(d_matrix_list), intent(in) :: dmat
    integer(IK), intent(in)         :: pw, depth
    integer(IK)                     :: i, j, bs
!
    this%ispc = 0
    this%uppd = 0
    do while (this%ispc < dmat%l)
      this%ispc = this%ispc + 1
      this%uppd = this%uppd + dmat%m(this%ispc)%g
      if (depth <= this%uppd) exit
    end do
    this%lowd = depth
    this%irow = dmat%m(this%ispc)%g - this%uppd + this%lowd
    this%nper = dmat%m(this%ispc)%g - this%irow + 1
    this%nsym = dmat%m(this%ispc)%s
    this%nnod = this%nper * this%nsym
    this%inod = 0
!
    allocate (this%nodes(this%nnod))
    bs = node__memsize(dmat)
!
    do concurrent(i=1:this%nper, j=1:this%nsym)
      block
        integer(IK) :: inod, ipnt
        inod = (j - 1) * this%nper + i
        ipnt = (inod - 1) * bs + pw
        call node__init(this%nodes(inod), i, j, ipnt)
      end block
    end do
!
    if (this%ispc < 1) return
!
  end subroutine breadth__init
!
  pure elemental function breadth__memsize(this, dmat) result(res)
    type(breadth_), intent(in)      :: this
    type(d_matrix_list), intent(in) :: dmat
    integer(IK)                     :: res
    res = node__memsize(dmat) * this%nnod
  end function breadth__memsize
!
  pure elemental subroutine breadth__destroy(this)
    type(breadth_), intent(inout) :: this
    if (ALLOCATED(this%nodes)) deallocate (this%nodes)
  end subroutine breadth__destroy
!
  pure subroutine tree_init(this, pw, blk, rot)
    class(tree), intent(inout)           :: this
    integer(IK), intent(in)              :: pw
    type(mol_block_list), intent(in)     :: blk
    type(molecular_rotation), intent(in), optional :: rot(*)
    integer(IK)                          :: n, i, ipw
!
    this%iscope = 0
!
    ipw = pw
    this%dmat = d_matrix_list(blk, ipw)
    ipw = ipw + this%dmat%memsize()
!
    call node__init_as_root(this%root, this%dmat, ipw)
    ipw = ipw + node__memsize(this%dmat)
!
    n = this%dmat%n_depth()
    if(n<1) return
!
    allocate (this%perms(n))
    allocate (this%breadthes(n))
!
    do i = 1, n
      call breadth__init(this%breadthes(i), this%dmat, ipw, i)
      ipw = ipw + breadth__memsize(this%breadthes(i), this%dmat)
    end do
!
    allocate (this%rots(this%dmat%l))
!
    if (PRESENT(rot)) then
      do concurrent(i=1:this%dmat%l)
        this%rots(i) = rot(i)
      end do
    end if
!
  end subroutine tree_init
!
  pure elemental function tree_memsize(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: res
    res = this%dmat%memsize() + node__memsize(this%dmat)
    if (.not. ALLOCATED(this%breadthes)) return
    res = res + SUM(breadth__memsize(this%breadthes, this%dmat))
  end function tree_memsize
!
  pure elemental function tree_n_depth(this) result(res)
    class(tree), intent(in) :: this
    integer(IK)             :: res
    if (ALLOCATED(this%breadthes)) then
      res = SIZE(this%breadthes)
    else
      res = 0
    end if
  end function tree_n_depth
!
   subroutine tree_setup(this, X, Y, W)
    class(tree), intent(inout) :: this
    real(RK), intent(in)       :: X(*)
    real(RK), intent(in)       :: Y(*)
    real(RK), intent(inout)    :: W(*)
    integer(IK)                :: nd
!
    if (.not. ALLOCATED(this%breadthes)) return
!
    call this%dmat%eval(this%rots, X, Y, W)
!
!!! lower bound is set to fixed lowerbound + Hungarian lowerbound.
!!! Assign the covariance matrix obtained from fixed points to root covariance matrix.
!
    W(this%root%v) = W(this%dmat%v) + W(this%dmat%o)
    W(this%root%h) = W(this%dmat%h)
    call copy(this%dmat%dd, W(this%dmat%c), W(this%root%c))
!
    nd = this%n_depth()
!
    if (nd < 1) then
      this%root%alive = .false.
      this%iscope = 0
      return
    else
      this%root%alive = .true.
    end if
!
    call init_perms(this%dmat, this%perms)
    call breadth__eval(this%breadthes(1), this%dmat, this%perms, &
   &                   W(this%root%h), W(this%root%c), W, nd == 1)
    this%iscope = 1
!
  end subroutine tree_setup
!
  pure subroutine init_perms(dmat, perms)
    type(d_matrix_list), intent(in) :: dmat
    integer(IK), intent(inout)      :: perms(*)
    integer(IK)                     :: i, j, k
    k = 0
    do j = 1, dmat%l
      do i = 1, dmat%m(j)%g
        k = k + 1
        perms(k) = i
      end do
    end do
  end subroutine init_perms
!
  pure elemental subroutine tree_destroy(this)
    type(tree), intent(inout) :: this
    call this%dmat%clear()
    if (ALLOCATED(this%breadthes)) deallocate (this%breadthes)
    if (ALLOCATED(this%perms)) deallocate (this%perms)
    if (ALLOCATED(this%rots)) deallocate (this%rots)
  end subroutine tree_destroy
!
!
  pure  subroutine breadth__eval(this, dmat, perm, H, C, W, is_terminal)
    type(breadth_), intent(inout)   :: this
    type(d_matrix_list), intent(in) :: dmat
    integer(IK), intent(in)         :: perm(*)
    real(RK), intent(in)            :: H, C(*)
    real(RK), intent(inout)         :: W(*)
    logical, intent(in)             :: is_terminal
    integer(IK)                     :: nr, i
!
    nr = this%uppd - this%lowd
    do concurrent(i=0:this%nper - 1)
      block
        integer(IK) :: j, ip, ir(nr)
!
        ip = perm(this%lowd + i)
        do concurrent(j=1:i)
          ir(j) = perm(this%lowd + j - 1)
        end do
        do concurrent(j=i + 1:nr)
          ir(j) = perm(this%lowd + j)
        end do
!
        do concurrent(j=1:this%nsym)
          block
            integer(IK) :: inod
!
            inod = this%nsym * i + j
            this%nodes(inod)%alive = .true.
            W(this%nodes(inod)%h) = H
            call copy(dmat%dd, C, W(this%nodes(inod)%c))
!
            call d_matrix_partial_eval(dmat%m(this%ispc), this%irow, ip, j, ir, W, &
           &                           W(this%nodes(inod)%v), W(this%nodes(inod)%h),  &
           &                           W(this%nodes(inod)%c))
!
          end block
        end do
      end block
    end do
!
    call breadth__set_node_index(this, W)
!
    if (is_terminal) then
      do concurrent(i=1:this%nnod)
        this%nodes(i)%alive = .false.
      end do
    end if
!
  end subroutine breadth__eval
!
  pure subroutine breadth__set_node_index(this, W)
    type(breadth_), intent(inout) :: this
    real(RK), intent(in)          :: W(*)
    real(RK)                      :: lv
    integer(IK)                   :: i
!
    this%inod = 0
!
    lv = RHUGE
    do i = 1, this%nnod
      if (this%nodes(i)%alive .and. W(this%nodes(i)%v) < lv) then
        this%inod = i
        lv = W(this%nodes(i)%v)
      end if
    end do
!
  end subroutine breadth__set_node_index
!
!
!!!
!
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
        real(RK)    :: H, C(1)
        integer(IK) :: ires(dmat%m(i)%g)
        integer(IK) :: j
        do concurrent(j=1:dmat%m(i)%g)
          ires(j) = j
        end do
        call d_matrix_partial_eval(dmat%m(i), 0, 0, 0, ires, W, res%L(i), H, C)
      end block
    end do
!
    res%lowerbound = W(dmat%v) + SUM(res%L)
!
    call next_index(dmat%l, dmat%m, 0, 0, res%ispc, res%depth)
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
    type(node)                      :: res
    integer(IK)                     :: i, np
!
    allocate (res%C, source=parent%C)
    res%H = parent%H
!
    np = SIZE(parent%per) - parent%depth
!
    block
      integer(IK) :: ip, ir(np)
!
      call perm_index(np, iper, parent%per(parent%depth), ip, ir)
      call d_matrix_partial_eval(dmat%m(parent%ispc), parent%depth, &
     &                           ip, isym, ir, W, res%lowerbound, res%H, res%C)
!
      res%lowerbound = res%lowerbound + SUM(parent%L)
!
      call next_index(dmat%l, dmat%m, parent%depth, parent%ispc, res%ispc, res%depth)
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
  pure subroutine next_index(l, m, depth, ispc, next_ispc, next_depth)
    integer(IK), intent(in)    :: l, depth, ispc
    type(d_matrix), intent(in) :: m(l)
    integer(IK), intent(inout) :: next_ispc, next_depth
    integer(IK)                :: i, spc
    spc = MAX(ispc, 1)
    if (depth < m(spc)%g) then
      next_ispc = spc
      next_depth = depth + 1
      return
    end if
    do i = spc + 1, l
      if (m(i)%g < 1) cycle
      next_ispc = i
      next_depth = 1
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
!!!
!
!| generate childe nodes instance
  pure function node_generate_breadth(this, dmat, W) result(res)
    class(node), intent(in)         :: this
    type(d_matrix_list), intent(in) :: dmat
    real(RK), intent(in)            :: W(*)
    type(breadth)                   :: res
    integer(IK)                     :: iper, isym
!
    if (.not. this%has_child()) then
      allocate (res%nodes(0))
      return
    end if
!
    res%nper = this%n_per(dmat)
    res%nsym = this%n_sym(dmat)
    res%nnod = res%nper * res%nsym
!
    allocate (res%nodes(res%nnod))
!
    do concurrent(iper=1:res%nper, isym=1:res%nsym)
      block
        integer(IK) :: inod
        inod = (isym - 1) * res%nper + iper
        res%nodes(inod) = node_new(this, iper, isym, dmat, W)
      end block
    end do
!
    call breadth_set_node_index(res)
!
  end function node_generate_breadth
!
  pure function breadth_generate_breadth(this, dmat, W) result(res)
    class(breadth), intent(in)      :: this
    type(d_matrix_list), intent(in) :: dmat
    real(RK), intent(in)            :: W(*)
    type(breadth)                   :: res
!
    if(this%inod<1)then
      allocate (res%nodes(0))
    else
      res = node_generate_breadth(this%nodes(this%inod), dmat, W)
    end if
!
  end function breadth_generate_breadth
!
  pure elemental subroutine breadth_prune(this, upperbound)
    class(breadth), intent(inout) :: this
    real(RK), intent(in)          :: upperbound
    integer(IK)                   :: i, n
    n = this%nper * this%nsym
    do concurrent(i=1:n)
      if (i == this%inod) then
        this%nodes(i)%alive = .false.
      else
        this%nodes(i)%alive = this%nodes(i)%alive .and. (this%nodes(i)%lowerbound < upperbound)
      end if
    end do
  end subroutine breadth_prune
!
  pure elemental function breadth_not_finished(this) result(res)
    class(breadth), intent(in) :: this
    logical                    :: res
    if (ALLOCATED(this%nodes)) then
      res = ANY(this%nodes%alive)
    else
      res = .false.
    end if
  end function breadth_not_finished
!
  pure elemental function breadth_is_finished(this) result(res)
    class(breadth), intent(in) :: this
    logical                    :: res
    res = .not. this%not_finished()
  end function breadth_is_finished
!
  pure elemental function breadth_iper(this) result(res)
    class(breadth), intent(in) :: this
    integer(IK)                :: res
    if (this%nnod < 1 .or. this%inod < 1) then
      res = 0
    else
      res = MODULO(this%inod - 1, this%nper)
    end if
  end function breadth_iper
!
  pure elemental function breadth_isym(this) result(res)
    class(breadth), intent(in) :: this
    integer(IK)                :: res
    if (this%nnod < 1 .or. this%inod < 1) then
      res = 0
    else
      res = 1 + (this%inod - 1) / this%nper
    end if
  end function breadth_isym
!
  pure elemental subroutine breadth_set_node_index(this, inod)
    class(breadth), intent(inout)     :: this
    integer(IK), intent(in), optional :: inod
!
    this%inod = 0
    if (.not.ALLOCATED(this%nodes)) return
!
    if (PRESENT(inod)) then
      if (0 < inod .and. inod <= this%nnod) then
        if(this%nodes(inod)%alive) this%inod = inod
        return
      end if
    end if
!
    block
      real(RK)    :: lv
      integer(IK) :: i
      lv = RHUGE
      do i = 1, this%nnod
        if (this%nodes(i)%alive .and. this%nodes(i)%lowerbound < lv) then
          this%inod = i
          lv = this%nodes(i)%lowerbound
        end if
      end do
    end block
!
  end subroutine breadth_set_node_index
!
  pure elemental subroutine breadth_set_node_minloc(this)
    class(breadth), intent(inout)     :: this
    if (.not.ALLOCATED(this%nodes)) return
    this%inod = MINLOC(this%nodes%lowerbound, 1)
  end subroutine breadth_set_node_minloc
!
  pure elemental function breadth_lowerbound(this) result(res)
    class(breadth), intent(in) :: this
    real(RK)                   :: res
    if (.not. ALLOCATED(this%nodes) .or. this%inod < 1) then
      res = RHUGE
    else
      res = this%nodes(this%inod)%lowerbound
    endif
  end function breadth_lowerbound
!
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
!!!
!
  pure subroutine copy(d, source, dest)
    integer(IK), intent(in) :: d
    real(RK), intent(in)    :: source(*)
    real(RK), intent(inout) :: dest(*)
    integer(IK)             :: i
    do concurrent(i=1:d)
      dest(i) = source(i)
    end do
  end subroutine copy
!
end module mod_tree
