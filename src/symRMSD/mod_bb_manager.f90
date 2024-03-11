!| Module for manage lowerbound matrices and tree.
!    By default, memory is allocated as follows.<br>
!    |-C1-|----W1----|<br>
!    |----|--C2--|--------W2-------|<br>
!    |-----------|--C3--|-W3-|<br>
!    Therefore, the maximum memory allocation size is MAX( SUM_i^I |Ci| + |W_I| ).
module mod_bb_manager
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_params, only: gemm, dot, copy
  use mod_mol_block
  use mod_c_matrix
  use mod_f_matrix
  use mod_rotation_matrix
  use mod_Hungarian
  use mod_tree
  implicit none
  private
  public :: bb_manager
  public :: bb_manager_memsize
  public :: bb_manager_worksize
  public :: bb_manager_setup
  public :: bb_manager_set_root
  public :: bb_manager_expand
  public :: bb_manager_select_top_node
  public :: bb_manager_leave
  public :: bb_manager_queue_is_empty
  public :: bb_manager_queue_is_bottom
  public :: bb_manager_current_value
!
  integer(IK), parameter :: header_size = 9
!
  integer(IK), parameter :: bq = 1
  !! pointer to mol_block interger array
  integer(IK), parameter :: cq = 2
  !! pointer to c_matrix interger array
  integer(IK), parameter :: fq = 3
  !! pointer to f_matrix interger array
  integer(IK), parameter :: tq = 4
  !! pointer to tree interger array
  integer(IK), parameter :: cx = 5
  !! pointer to c_matrix memory
  integer(IK), parameter :: fx = 6
  !! pointer to f_matrix memory
  integer(IK), parameter :: tx = 7
  !! pointer to tree memory
  integer(IK), parameter :: cw = 8
  !! pointer to c_matrix work memory
  integer(IK), parameter :: fw = 9
  !! pointer to f_matrix work memory
!
!| bb_manager<br>
!    Node data [L, G, C, F, W]<br>
!    0    L      : scalar, lowerbound.<br>
!    1    G      : scalar, partial sum of auto variance.<br>
!    2    F(n,n) : Free rotate score matrix.<br>
!    2+nn C(D,D) : partial sum of covariance.<br>
!         W(*)   : work space.<br>
!  This is mainly used for passing during initialization.
  type bb_manager
    integer(IK), allocatable :: q(:)
    !! work integer array
    real(RK), allocatable    :: x(:)
    !! main memory
    real(RK), allocatable    :: w(:)
    !! work array
  contains
    final           :: bb_manager_destroy
  end type bb_manager
!
  interface bb_manager
    module procedure bb_manager_new
  end interface bb_manager
!
contains
!
!| Constructer
  pure function bb_manager_new(m, n, sym) result(res)
    integer(IK), intent(in) :: m
    !! number of molecules.
    integer(IK), intent(in) :: n
    !! number of atoms per molecule.
    integer(IK), intent(in), optional :: sym(:,:)
    !! symmetric codomains, [[a1,a2,...,am], [b1,b2,...,bm], ...].
    type(mol_block)         :: b
    type(c_matrix)          :: c
    type(f_matrix)          :: f
    type(tree)              :: t
    type(bb_manager)        :: res
    integer(IK)             :: q(header_size)
!
    b = mol_block(m, n, sym)
    c = c_matrix(b%q)
    f = f_matrix(b%q)
    t = tree(b%q, node_memsize)
!
    q(bq) = header_size + 1
    q(cq) = q(bq) + SIZE(b%q)
    q(fq) = q(cq) + SIZE(c%q)
    q(tq) = q(fq) + SIZE(f%q)
!
    q(cx) = 1
    q(cw) = q(cx) + c_matrix_memsize(c%q)
    q(tx) = q(cw)
    q(fx) = q(cw) + 2 + DD                 ! F as root node
    q(fw) = q(fx) + f_matrix_memsize(f%q)
!
    allocate (res%q, source=[q, b%q, c%q, f%q, t%q])
    allocate (res%x(bb_manager_memsize(res%q)))
    allocate (res%w(bb_manager_worksize(res%q)))
!
  end function bb_manager_new
!
  pure function node_memsize(b, p) result(res)
    integer(IK), intent(in) :: b(*)
    integer(IK), intent(in) :: p
    integer(IK)             :: res, n
    n = mol_block_nmol(b) - p
    res = 1 + 1 + DD + n * n
  end function node_memsize
!
  pure function tree_worksize(b) result(res)
    integer(IK), intent(in) :: b(*)
    integer(IK)             :: p, n, m, sw, hw
    integer(IK)             :: res
    sw = sdmin_worksize()
    res = 1 + mol_block_nsym(b) * sw
    n = mol_block_nmol(b)
    m = n
    do p = 1, n
      hw = Hungarian_worksize(m, m)
      res = MAX(hw, res)
      m = m - 1
    end do
  end function tree_worksize
!
!| Inquire worksize of f_matrix.
  pure function bb_manager_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_manager.
    integer(IK)             :: res
    res = c_matrix_memsize(q(q(cq))) + tree_memsize(q(q(tq)))
  end function bb_manager_memsize
!
!| Inquire worksize of bb_manager.
  pure function bb_manager_worksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array.
    integer(IK)             :: res
    res = MAX(tree_worksize(q(q(bq))), &
   &          c_matrix_worksize(q(q(cq))) - tree_memsize(q(q(tq))) - tree_worksize(q(q(bq))))
  end function bb_manager_worksize
!
!| Setup C matrix and F matrix in root node.
  pure subroutine bb_manager_setup(Q, X, Y, W)
    integer(IK), intent(in) :: Q(*)
    !! integer array
    real(RK), intent(in)    :: X(*)
    !! reference coordinate
    real(RK), intent(in)    :: Y(*)
    !! target coordinate
    real(RK), intent(inout) :: W(*)
    !! work integer array
!
    call c_matrix_eval(Q(Q(cq)), Q(Q(bq)), X, Y, W(Q(cx)), W(Q(cw)))
    call f_matrix_eval(Q(Q(fq)), Q(Q(cq)), W(Q(cx)), W(Q(fx)), W(Q(fw)))
    call zfill(DD + 2, W(Q(tx)), 1)
!
  end subroutine bb_manager_setup
!
!| Setup C matrix and F matrix in root node.
  pure subroutine bb_manager_set_root(QP, XP, Q, X)
    integer(IK), intent(in)      :: QP(*)
    !! parent integer array
    real(RK), intent(in)         :: XP(*)
    !! parent work array
    integer(IK), intent(in)      :: Q(*)
    !! integer array
    real(RK), intent(inout)      :: X(*)
    !! work array
    integer(IK)                  :: pc, pp
!
    pc = Q(Q(tx)) + 1
    pp = tree_current_pointer(QP(QP(tq))) + 1
    X(pc) = XP(pp)
    pc = pc + 1
    pp = pp + 1
    call copy(DD, XP(pp), 1, X(pc), 1)
!
  end subroutine bb_manager_set_root
!
!| Expand top node in queue.
  subroutine bb_manager_expand(Q, X, W)
    integer(IK), intent(inout)   :: Q(*)
    !! work integer array
    real(RK), intent(inout)      :: X(*)
    !! main memory
    real(RK), intent(out)        :: W(*)
    !! work real array
    integer(IK)                  :: nper, nsym
    integer(IK)                  :: p, n, np, nn, nb, nw
!
     np = Q(tx) - 1 + tree_current_pointer(Q(Q(tq)))
     call tree_expand(Q(Q(tq)))
!
     p = tree_current_level(Q(Q(tq)))
     n = mol_block_nmol(Q(Q(bq)))
     nn = Q(tx) - 1 + tree_queue_pointer(Q(Q(tq)))
     nb = node_memsize(Q(Q(bq)), p)
     nper = tree_n_perm(Q(Q(tq)))
     nsym = mol_block_nsym(Q(Q(bq)))
     nw = MAX(Hungarian_worksize(n - p, n - p), sdmin_worksize())
!
     block
       integer(IK) :: s(n)
       s = tree_current_permutation(Q(Q(tq)))
       print*,s
       if (n == p) then
         call expand_terminal(p, n, nb, nw, nsym, Q(Q(cq)), Q(Q(bq)), s, X(Q(cx)), X(np), X(nn), W)
       else
         call expand(p, n, nb, nw, nper, nsym, Q(Q(cq)), Q(Q(bq)), s, X(Q(cx)), X(np), X(nn), W(1), W(2))
       end if
     end block
!
  end subroutine bb_manager_expand
!
  pure subroutine expand(p, n, nb, nw, nper, nsym, cq, b, s, C, NP, NN, W1, W2)
    integer(IK), intent(in)     :: p, n, nb, nw, nper, nsym
    integer(IK), intent(in)     :: cq(*), b(*), s(*)
    real(RK), intent(in)        :: C(*)
    real(RK), intent(in)        :: NP(*)
    real(RK), intent(inout)     :: NN(nb, nsym, nper)
    real(RK), intent(out)       :: W1(*)
    real(RK), intent(out)       :: W2(nw, *)
    integer(IK), parameter      :: mmap_L = 1
    integer(IK), parameter      :: mmap_G = 2
    integer(IK), parameter      :: mmap_C = 3
    integer(IK)                 :: mmap_F
    integer(IK)                 :: mp, mn, mm
    integer(IK)                 :: iper, isym
!
    mn = n - p
    mm = mn**2
    mp = mn + 1
!
    mmap_F = mmap_C + DD
!
    do concurrent(iper=1:nper, isym=1:nsym)
      call copy(DD + 1, NP(mmap_G), 1, NN(mmap_G, isym, iper), 1)
      call c_matrix_add(cq, p, s(iper + p - 1), isym, C, NN(mmap_G, isym, iper), NN(mmap_C, isym, iper))
    end do
!
    do iper = 1, nper
!
        call subm(mp, mp, iper, NP(mmap_F), NN(mmap_F, 1, iper))
        call Hungarian(mn, mn, NN(mmap_F, 1, iper), w1)
        NN(mmap_L, 1, iper) = w1(1)
        call estimate_sdmin(NN(mmap_G, 1, iper), NN(mmap_C, 1, iper), w1)
!
        do concurrent(isym=2:nsym)
            call copy(mm, NN(mmap_F, 1, iper), 1, NN(mmap_F, isym, iper), 1)
            call estimate_sdmin(NN(mmap_G, isym, iper), NN(mmap_C, isym, iper), W2(1, isym - 1))
            NN(mmap_L, isym, iper) = W2(1, isym - 1) + NN(mmap_L, 1, iper)
        end do
!
        NN(mmap_L, 1, iper) = W1(1) + NN(mmap_L, 1, iper)
!
    end do
!
  end subroutine expand
!
  pure subroutine expand_terminal(p, n, nb, nw, nsym, cq, b, s, C, NP, NN, W)
    integer(IK), intent(in)     :: p, n, nb, nw, nsym
    integer(IK), intent(in)     :: cq(*), b(*), s(*)
    real(RK), intent(in)        :: C(*)
    real(RK), intent(in)        :: NP(*)
    real(RK), intent(inout)     :: NN(nb, nsym)
    real(RK), intent(out)       :: W(nw, *)
    integer(IK), parameter      :: mmap_L = 1
    integer(IK), parameter      :: mmap_G = 2
    integer(IK), parameter      :: mmap_C = 3
    integer(IK)                 :: isym
!
    do concurrent(isym=1:nsym)
      call copy(DD + 1, NP(mmap_G), 1, NN(mmap_G, isym), 1)
      call c_matrix_add(cq, p, s(n), isym, C, NN(mmap_G, isym), NN(mmap_C, isym))
      call estimate_sdmin(NN(mmap_G, isym), NN(mmap_C, isym), W(1, isym))
      NN(mmap_L, isym) = W(1, isym)
    end do
!
  end subroutine expand_terminal
!
  pure subroutine subm(m, n, r, X, Y)
  integer(IK), intent(in) :: m, n, r
  real(RK), intent(in)    :: X(m, n)
  real(RK), intent(inout) :: Y(m - 1, n - 1)
  integer(IK)             :: i, j
    do concurrent(j=2:n)
      do concurrent(i=1:r - 1)
        Y(i, j - 1) = X(i, j)
      end do
      do concurrent(i=r + 1:m)
        Y(i - 1, j - 1) = X(i, j)
      end do
    end do
  end subroutine subm
!
!| Select_top_node.
  subroutine bb_manager_select_top_node(Q, X, UB)
    integer(IK), intent(inout) :: Q(*)
    !! work integer array
    real(RK), intent(in)       :: X(*)
    !! main memory
    real(RK), intent(in)       :: UB
    !! upper bound
    call tree_select_top_node(Q(Q(tq)), UB, X(Q(tx)))
     print*,'select top node',Q(Q(tq):Q(tq)+3)
     print'(4i4)',Q(Q(tq)+4:Q(tq)+19)
  end subroutine bb_manager_select_top_node
!
!| Leave current node.
  pure function bb_manager_queue_is_empty(Q) result(res)
    integer(IK), intent(in) :: Q(*)
    !! work integer array
    logical                 :: res
    res = tree_queue_is_empty(Q(Q(tq)))
  end function bb_manager_queue_is_empty
!
!| Leave current node.
  pure function bb_manager_current_value(Q, X) result(res)
    integer(IK), intent(in) :: Q(*)
    !! work integer array
    real(RK), intent(in)    :: X(*)
    !! main memory
    real(RK)                :: res
    integer(IK)             :: t
    t = Q(tx) - 1 + tree_current_pointer(Q(Q(tq)))
    res = X(t)
  end function bb_manager_current_value
!
!| Leave current node.
  pure function bb_manager_queue_is_bottom(Q) result(res)
    integer(IK), intent(in) :: Q(*)
    !! work integer array
    logical                 :: res
    res = tree_queue_is_bottom(Q(Q(tq)))
  end function bb_manager_queue_is_bottom
!
!| Leave current node.
  pure subroutine bb_manager_leave(Q)
    integer(IK), intent(inout)   :: Q(*)
    !! work integer array
!
    call tree_leave(Q(Q(tq)))
!
  end subroutine bb_manager_leave
!
! pure subroutine eval_child(b, p, inode, Wp, Wc)
!   type(mol_block), intent(in) :: b
!   integer(IK), intent(in)     :: p, inode
!   real(RK), intent(in)        :: Wp(*)
!   real(RK), intent(inout)     :: Wc(*)
!   integer(IK), parameter      :: l = 1
!   integer(IK), parameter      :: f = 2
!   integer(IK), parameter      :: g = 3
!   integer(IK)                 :: c
!
!   c = g + b%n1 * b%n2
!   call copy(1 + DD, Wp(c), 1, Wc(c), 1)
!   call subm(b%n1, b%n2, inode, Wp(f), Wc(f))
!
! contains
!
!   pure subroutine subm(m, n, r, X, Y)
!   integer(IK), intent(in) :: m, n, r
!   real(RK), intent(in)    :: X(m, n)
!   real(RK), intent(inout) :: Y(m-1, n-1)
!   integer(IK)             :: i, j
!     do concurrent(j = 2: n)
!       do concurrent(i=1:r - 1)
!         Y(i, j - 1) = X(i, j)
!       end do
!       do concurrent(i=r + 1:m)
!         Y(i - 1, j - 1) = X(i, j)
!       end do
!     end do
!   end subroutine subm

! end subroutine eval_child
!
! pure subroutine d_matrix_partial_eval(a, p, iprm, isym, ires, W, LT, H, C, LF, LB)
!   type(d_matrix), intent(in)        :: a
!   integer(IK), intent(in)           :: p, iprm, isym, ires(*)
!   real(RK), intent(in)              :: W(*)
!   real(RK), intent(inout)           :: LT, H, C(*)
!   real(RK), intent(inout), optional :: LF, LB
!
!   LT = ZERO
!   if (PRESENT(LF)) LF = ZERO
!   if (PRESENT(LB)) LB = ZERO
!
!   if (p < 0 .or. a%g < p) return
!   if (p < a%g) call setminus_eval(a, p, ires, W, LT)
!   if (PRESENT(LB)) LB = LT
!
!   if (p == 0) return
!   if (iprm < 1 .or. a%g < iprm) return
!   if (isym < 0 .or. a%s <= isym) return
!
!   block
!     integer(IK) :: ih, ic
!     ih = a%c + (iprm - 1) * a%cb + (p - 1) * a%cl
!     ic = ih + 1 + DD * isym
!     if (PRESENT(LF)) then
!       call partial_eval(a%s, a%g, a%nw2, p, W(ih), W(ic), LF, H, C)
!       LT = LT + LF
!     else
!       call partial_eval(a%s, a%g, a%nw2, p, W(ih), W(ic), LT, H, C)
!     end if
!   end block
!
! contains
!
!   pure subroutine setminus_eval(a, p, ires, W, res)
!     type(d_matrix), intent(in) :: a
!     integer(IK), intent(in)    :: p, ires(*)
!     real(RK), intent(in)       :: W(*)
!     real(RK), intent(inout)    :: res
!     integer(IK)                :: nw, nr, rr
!
!     if (p < 0 .or. a%g < p) then
!       res = ZERO
!     elseif (p == 0) then
!       nw = nwork(a%g)
!       block
!         real(RK) :: T(nw)
!         call Hungarian(a%g, a%g, W(a%z), T(1))
!         res = T(1)
!       end block
!     else
!       nr = a%g - p
!       rr = nr * nr
!       nw = rr + nwork(nr)
!       block
!         integer(IK) :: iw, iz
!         real(RK)    :: T(nw)
!         iw = 1 + rr
!         iz = a%z + a%g * (p - 1)
!         call pack_Z(a%g, nr, ires, W(a%z), T(1))
!         call Hungarian(nr, nr, T(1), T(iw))
!         res = T(iw)
!       end block
!     end if
!
!   end subroutine setminus_eval
!
!   pure subroutine partial_eval(s, n, nw, p, H, C, LF, HP, CP)
!     integer(IK), intent(in) :: s, n, nw
!     integer(IK), intent(in) :: p
!     real(RK), intent(in)    :: H, C(*)
!     real(RK), intent(inout) :: LF, HP, CP(*)
!     real(RK)                :: W(nw)
!     integer(IK), parameter  :: it = 1
!     integer(IK), parameter  :: ic = 2
!     integer(IK)             :: iw
!
!     iw = ic + DD
!
!!! update H and C
!     w(it) = HP + H
!     HP = w(it)
!     call add(DD, CP, C, W(ic))
!     call DCOPY(DD, W(ic), 1, CP, 1)
!
!!! get squared displacement
!     call estimate_sdmin(w(it), w(ic), W(iw))
!     LF = LF + w(iw)
!
!   end subroutine partial_eval
!
!   pure elemental function nwork(np) result(res)
!     integer(IK), intent(in) :: np
!     real(RK)                :: W(1)
!     integer(IK)             :: res
!     call Hungarian(-np, -np, W, W)
!     res = NINT(W(1), IK)
!   end function nwork
!
!   pure subroutine pack_Z(n, nr, ires, Z, T)
!     integer(IK), intent(in) :: n, nr, ires(nr)
!     real(RK), intent(in)    :: Z(n, *)
!     real(RK), intent(inout) :: T(nr, nr)
!     integer(IK)             :: i, j
!     do concurrent(i=1:nr, j=1:nr)
!       T(i, j) = Z(ires(i), n - nr + j)
!     end do
!   end subroutine pack_Z
!
! end subroutine d_matrix_partial_eval
!
!!! d_matrix_list
!
!| Constructor of d_matrix_list
! pure function d_matrix_list_new(b, p) result(res)
!   type(mol_block_list), intent(in)     :: b
!   integer(IK), intent(in)              :: p
!   type(d_matrix_list)                  :: res
!   integer(IK)                          :: i, ip
!
!   res%l = b%nspecies()
!   res%h = p              ! H(1)
!   res%v = res%h + 1      ! V(1)
!   res%c = res%v + 1      ! C(d*d)
!   res%o = res%c + DD     ! O(L+1)
!   ip = res%o + res%l + 1
!
!   res%nk = worksize_sdmin()
!
!   allocate (res%m(res%l))
!
!   do i = 1, res%l
!     res%m(i) = d_matrix(ip, res%nk, b%b(i))
!     ip = ip + d_matrix_memsize(res%m(i))
!   end do
!
! end function d_matrix_list_new
!
! pure function d_matrix_list_memsize(this) result(res)
!   class(d_matrix_list), intent(in) :: this
!   integer(IK)                      :: res
!   if (ALLOCATED(this%m)) then
!     res = SUM(d_matrix_memsize(this%m)) + (this%l + 1) + DD + 2
!   else
!     res = 0
!   end if
! end function d_matrix_list_memsize
!
! pure elemental function d_matrix_list_n_depth(this) result(res)
!   class(d_matrix_list), intent(in) :: this
!   integer(IK)                      :: res
!   if (ALLOCATED(this%m)) then
!     res = SUM(this%m%g)
!   else
!     res = 0
!   end if
! end function d_matrix_list_n_depth
!
! pure subroutine d_matrix_list_eval(this, rot, X, Y, W)
!   class(d_matrix_list), intent(in) :: this
!   type(mol_symmetry), intent(in)   :: rot(*)
!   real(RK), intent(in)             :: X(*), Y(*)
!   real(RK), intent(inout)          :: W(*)
!   integer(IK)                      :: i
!
!   if (.not. ALLOCATED(this%m)) return
!
!!! estimate H_fix, V_fix, C_fix. (mol_blocks for g<n or (g=1, s=1))
!!! if set is empty, H_fix=0, V_fix=0, C_fix=0.
!
!   call fixpoints_eval(this%nk, this%l, this%m, X, Y,  &
!  &                    W(this%h), W(this%v), W(this%c))
!
!!! estimate floating mol blocks.
!!! O is filled by lower bounds for each blocks.
!
!   do concurrent(i=1:this%l)
!     call d_matrix_eval(this%m(i), rot(i), X, Y, W)
!     block
!       integer(IK) :: j, ires(this%m(i)%g)
!       do concurrent(j=1:this%m(i)%g)
!         ires(j) = j
!       end do
!       j = this%o + i - 1
!       call d_matrix_partial_eval(this%m(i), 0, 0, 0, ires, W, W(j), W(j), W(j))
!     end block
!   end do
!
!!! transform O to accumulated form.
!
!   do concurrent(i=this%l - 1:1:-1)
!     W(this%o + i - 1) = W(this%o + i - 1) + W(this%o + i)
!   end do
!
!   W(this%o + this%l) = ZERO
!
! contains
!
!   pure subroutine fixpoints_eval(nk, l, m, X, Y, H, V, C)
!     integer(IK), intent(in)    :: nk, l
!     type(d_matrix), intent(in) :: m(l)
!     real(RK), intent(in)       :: X(*), Y(*)
!     real(RK), intent(inout)    :: H, V, C(*)
!     integer(IK)                :: t(l), p(l), q(l)
!     integer(IK), parameter     :: ig = 1
!     integer(IK), parameter     :: iv = 2
!     integer(IK), parameter     :: ic = 3
!     integer(IK)                :: i, iw, ix, iy, mn, dmn, dmn2, nw
!
!     do concurrent(i=1:l)
!       t(i) = m(i)%m * (m(i)%n - m(i)%g)
!     end do
!
!     mn = SUM(t)
!     dmn = D * mn
!     dmn2 = dmn + dmn
!
!     iw = ic + DD
!     ix = ic + DD
!     iy = ix + dmn
!     nw = 2 + DD + MAX(dmn + dmn, DD + nk)
!
!     do concurrent(i=1:l)
!       t(i) = t(i) * D
!     end do
!     p(1) = 0
!     do i = 2, l
!       p(i) = p(i - 1) + t(i - 1)
!     end do
!     do concurrent(i=1:l)
!       q(i) = m(i)%x + m(i)%dm * m(i)%g
!     end do
!
!     block
!       real(RK) :: W(nw)
!
!       do concurrent(i=1:l)
!         block
!           integer(IK) :: px
!           px = p(i) + ix
!           call DCOPY(t(i), X(q(i)), 1, W(px), 1)
!         end block
!       end do
!       do concurrent(i=1:l)
!         block
!           integer(IK) :: py
!           py = p(i) + iy
!           call DCOPY(t(i), Y(q(i)), 1, W(py), 1)
!         end block
!       end do
!
!       W(ig) = ddot(dmn2, W(ix), 1, W(ix), 1)
!
!       if (mn > 0) then
!         call DGEMM('N', 'T', D, D, mn, ONE, W(iy), D, W(ix), D, ZERO, W(ic), D)
!         call estimate_sdmin(W(ig), W(ic), W(iw))
!         w(iv) = w(iw)
!       else
!         call zfill(DD, W(ic))
!         w(iv) = ZERO
!       end if
!
!       H = W(ig)
!       V = W(ig) - W(iv) - W(iv)
!       call DCOPY(DD, W(ic), 1, C, 1)
!
!     end block
!
!   end subroutine fixpoints_eval
!
! end subroutine d_matrix_list_eval
!
! pure subroutine d_matrix_list_partial_eval(this, p, perm, iprm, isym, W, LT, H, C, LF, LB)
!   class(d_matrix_list), intent(in)  :: this
!   integer(IK), intent(in)           :: p, perm(*), iprm, isym
!   real(RK), intent(in)              :: W(*)
!   real(RK), intent(inout)           :: LT, H, C(*)
!   real(RK), intent(inout), optional :: LF, LB
!   integer(IK)                       :: ispc, iofs, nres
!
!   call p_index(this, p, ispc, iofs)
!   if (ispc < 1) return
!
!   nres = this%m(ispc)%g - iofs
!   block
!     integer(IK) :: i, jper, ires(nres)
!     jper = perm(p - 1 + iprm)
!     do concurrent(i=1:iprm - 1)
!       ires(i) = perm(p + i - 1)
!     end do
!     do concurrent(i=iprm:nres)
!       ires(i) = perm(p + i)
!     end do
!     call d_matrix_partial_eval(this%m(ispc), iofs, jper, isym, ires, W, LT, H, C, LF, LB)
!   end block
!
! contains
!
!   pure elemental subroutine p_index(this, p, ispc, iofs)
!     type(d_matrix_list), intent(in) :: this
!     integer(IK), intent(in)         :: p
!     integer(IK), intent(inout)      :: ispc, iofs
!     integer(IK)                     :: i, q, r
!
!     ispc = 0
!     iofs = 0
!     r = p
!     q = 1
!
!     do i = 1, this%l
!       if (q > p) return
!       ispc = i
!       iofs = r
!       r = r - this%m(i)%g
!       q = q + this%m(i)%g
!     end do
!
!     if (q <= p) ispc = 0
!
!   end subroutine p_index
!
! end subroutine d_matrix_list_partial_eval
!
! pure elemental subroutine d_matrix_list_clear(this)
!   class(d_matrix_list), intent(inout) :: this
!   this%l = 0
!   if (ALLOCATED(this%m)) deallocate (this%m)
! end subroutine d_matrix_list_clear
!
! pure elemental subroutine d_matrix_list_destroy(this)
!   type(d_matrix_list), intent(inout) :: this
!   call d_matrix_list_clear(this)
! end subroutine d_matrix_list_destroy
!
!!! util
!
! pure subroutine add(d, A, B, C)
!   integer(IK), intent(in) :: d
!   real(RK), intent(in)    :: A(*), B(*)
!   real(RK), intent(inout) :: C(*)
!   integer(IK)             :: i
!   do concurrent(i=1:d)
!     C(i) = A(i) + B(i)
!   end do
! end subroutine add
!
  pure subroutine zfill(d, x, ld)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: x(*)
    integer(IK), intent(in) :: ld
    integer(IK)             :: i, dld
    dld = d * ld
    do concurrent(i=1:dld:ld)
      x(i) = ZERO
    end do
  end subroutine zfill
!
  pure elemental subroutine bb_manager_destroy(this)
    type(bb_manager), intent(inout) :: this
    if (ALLOCATED(this%x)) deallocate (this%x)
    if (ALLOCATED(this%w)) deallocate (this%w)
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine bb_manager_destroy
!
!| Setup.
! pure subroutine bb_manager_list_setup(bb, X, Y, W)
!   type(bb_manager), intent(in) :: bb(:)
!   !! bb_manager
!   real(RK), intent(in)         :: X(*)
!   !! reference coordinate
!   real(RK), intent(in)         :: Y(*)
!   !! target coordinate
!   real(RK), intent(inout)      :: W(*)
!   !! work memory
!   integer(IK)                  :: k, l
!
!   do k = 1, SIZE(bb)
!     call c_matrix_eval(bb(k)%c, bb(k)%b, bb(k)%ms, X, Y, W, W)
!     call f_matrix_eval(bb(k)%f, bb(k)%b, W(bb(k)%c%p), W, W)
!     call Hungarian(bb(k)%b%n1, bb(k)%b%n2, W(bb(k)%f%p), W(bb(k)%t%p))
!     do concurrent(l=1:k - 1)
!       W(bb(l)%o) = W(bb(l)%o) + W(bb(k)%t%p)
!     end do
!     W(bb(k)%o) = ZERO
!     call copy(bb(k)%b%n1 * bb(k)%b%n2, W(bb(k)%f%p), 1, W(bb(k)%t%p + mmap_f), 1)
!     call zfill(DD + 1, W(bb(k)%t%p + mmap_g(bb(k)%b, 0)))
!   end do
!
! end subroutine bb_manager_list_setup
!
end module mod_bb_manager

