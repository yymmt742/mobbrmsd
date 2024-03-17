!| Module for manage lowerbound matrices and tree.
!    By default, memory is allocated as follows.<br>
!    |-C1-|----W1----|<br>
!    |----|--C2--|--------W2-------|<br>
!    |-----------|--C3--|-W3-|<br>
!    Therefore, the maximum memory allocation size is MAX( SUM_i^I |Ci| + |W_I| ).
module mod_bb_block
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
  public :: bb_block
  public :: bb_block_molsize
  public :: bb_block_memsize
  public :: bb_block_worksize
  public :: bb_block_setup
  public :: bb_block_set_root
  public :: bb_block_expand
  public :: bb_block_leave
  public :: bb_block_queue_is_empty
  public :: bb_block_queue_is_bottom
  public :: bb_block_current_value
!
  integer(IK), parameter :: header_size = 8
!
  integer(IK), parameter :: cx = 1
  !! pointer to c_matrix array (fixed)
  integer(IK), parameter :: ts = 1
  !! pointer to tree interger work array (fixed)
  integer(IK), parameter :: bq = header_size + 1
  !! pointer to mol_block interger array (fixed)
!
  integer(IK), parameter :: cq = 1
  !! pointer to c_matrix interger array
  integer(IK), parameter :: fq = 2
  !! pointer to f_matrix interger array
  integer(IK), parameter :: tq = 3
  !! pointer to tree interger array
  integer(IK), parameter :: ps = 4
  !! pointer to tree interger array
  integer(IK), parameter :: fx = 5
  !! pointer to f_matrix memory
  integer(IK), parameter :: tx = 6
  !! pointer to tree memory
  integer(IK), parameter :: cw = 7
  !! pointer to c_matrix work memory
  integer(IK), parameter :: fw = 8
  !! pointer to f_matrix work memory
!
!| bb_block<br>
!    Node data [L, G, C, F, W]<br>
!    0    L      : scalar, lowerbound.<br>
!    1    G      : scalar, partial sum of auto variance.<br>
!    2    F(n,n) : Free rotate score matrix.<br>
!    2+nn C(D,D) : partial sum of covariance.<br>
!         W(*)   : work space.<br>
!  This is mainly used for passing during initialization.
  type bb_block
    integer(IK), allocatable :: q(:)
    !! static integer array
    integer(IK), allocatable :: s(:)
    !! work integer array
  contains
    final           :: bb_block_destroy
  end type bb_block
!
  interface bb_block
    module procedure bb_block_new
  end interface bb_block
!
contains
!
!| Constructer
  pure function bb_block_new(m, n, sym) result(res)
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
    type(bb_block)          :: res
    integer(IK)             :: q(header_size)
!
    b = mol_block(m, n, sym)
    c = c_matrix(b%q)
    f = f_matrix(b%q)
    t = tree(b%q, node_memsize)
!
    q(cq) = bq + SIZE(b%q)
    q(fq) = q(cq) + SIZE(c%q)
    q(tq) = q(fq) + SIZE(f%q)
!
    q(ps) = ts + SIZE(t%s)
!
    q(cw) = cx + c_matrix_memsize(c%q)
    q(tx) = q(cw)
    q(fx) = q(cw) + 2 + DD ! F as root node
    q(fw) = q(fx) + f_matrix_memsize(f%q)
!
    allocate (res%q, source=[q, b%q, c%q, f%q, t%q])
    block
      integer(IK) :: prm(mol_block_nmol(b%q))
      call perm_reset(b%q, prm)
      allocate (res%s, source=[t%s, prm])
    end block
!
  end function bb_block_new
!
  pure function node_memsize(b, p) result(res)
    integer(IK), intent(in) :: b(*)
    integer(IK), intent(in) :: p
    integer(IK)             :: res, n
!
    n = mol_block_nmol(b) - p
    res = 1 + 1 + DD + n * n  ! [L, G, C(D,D), F(n,n)]
!
  end function node_memsize
!
  pure function tree_worksize(b) result(res)
    integer(IK), intent(in) :: b(*)
    integer(IK)             :: p, n, m, sw, hw
    integer(IK)             :: res
!
    sw = sdmin_worksize()
    res = 1 + mol_block_nsym(b) * sw
    n = mol_block_nmol(b)
    m = n
    do p = 1, n
      hw = Hungarian_worksize(m, m)
      res = MAX(hw, res)
      m = m - 1
    end do
!
  end function tree_worksize
!
!| Returns the memory size array size.
  pure function bb_block_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block header array.
    integer(IK)             :: res
!
    res = c_matrix_memsize(q(q(cq))) + tree_memsize(q(q(tq)))
!
  end function bb_block_memsize
!
!| Returns the memory size of work array size.
  pure function bb_block_worksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block header array.
    integer(IK)             :: res
!
    res = tree_worksize(q(bq))
    res = MAX(res, c_matrix_worksize(q(q(cq))) - tree_memsize(q(q(tq))) - res)
!
  end function bb_block_worksize
!
!| Inquire Returns the memory size of molecular block size.
  pure function node_molsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block header array.
    integer(IK)             :: res, n
    res = mol_block_total_size(q(bq))
  end function node_molsize
!
!| Setup C matrix and F matrix in root node.
  pure subroutine bb_block_setup(q, X, Y, s, W)
    integer(IK), intent(in)    :: q(*)
    !! integer array
    real(RK), intent(in)       :: X(*)
    !! reference coordinate
    real(RK), intent(in)       :: Y(*)
    !! target coordinate
    integer(IK), intent(inout) :: s(*)
    !! integer work array
    real(RK), intent(inout)    :: W(*)
    !! work integer array
!
    call tree_reset(q(q(tq)), s(ts))
    call perm_reset(q(bq), s(q(ps)))
    call c_matrix_eval(q(q(cq)), q(bq), X, Y, W(cx), W(q(cw)))
    call f_matrix_eval(q(q(fq)), q(q(cq)), W(cx), W(q(fx)), W(q(fw)))
    call zfill(DD + 2, W(q(tx)), 1)
!
  end subroutine bb_block_setup
!
  pure subroutine perm_reset(b, prm)
    integer(IK), intent(in)    :: b(*)
    !! mol_block
    integer(IK), intent(inout) :: prm(*)
    !! permutation
    integer(IK)                :: i, n
!
    n = mol_block_nmol(b)
    do concurrent(i=1:n)
      prm(i) = i
    end do
!
  end subroutine perm_reset
!
!| Setup C matrix and F matrix in root node.
  pure subroutine bb_block_set_root(qp, sp, XP, q, s, X)
    integer(IK), intent(in)      :: qp(*)
    !! parent integer array
    integer(ik), intent(in)      :: sp(*)
    !! parent integer work array
    real(RK), intent(in)         :: XP(*)
    !! parent work array
    integer(ik), intent(in)      :: q(*)
    !! integer array
    integer(ik), intent(inout)   :: s(*)
    !! integer work array
    real(RK), intent(inout)      :: X(*)
    !! work array
    integer(IK)                  :: pc, pp
!
    pc = q(q(tx)) + 1
    pp = tree_current_pointer(qp(qp(tq)), sp) + 1
    X(pc) = XP(pp)
    pc = pc + 1
    pp = pp + 1
    call copy(DD, XP(pp), 1, X(pc), 1)
!
    call tree_reset(q(q(tq)), s(ts))
    call perm_reset(q(bq), s(q(ps)))
!
  end subroutine bb_block_set_root
!
!| Expand top node in queue.
  pure subroutine bb_block_expand(q, UB, s, X, W)
    integer(IK), intent(in)    :: q(*)
    !! work integer array
    real(RK), intent(in)       :: UB
    !! upper bound
    integer(IK), intent(inout) :: s(*)
    !! work integer array
    real(RK), intent(inout)    :: X(*)
    !! main memory
    real(RK), intent(inout)    :: W(*)
    !! work real array
    integer(IK)                :: nper, nsym
    integer(IK)                :: p, n, np, nn, nb, nw
!
     np = q(tx) - 1 + tree_current_pointer(q(q(tq)), s(ts))
     call tree_expand(q(q(tq)), s(ts))
!
     p = tree_current_level(s)
     n = mol_block_nmol(q(bq))
     nn = q(tx) - 1 + tree_queue_pointer(q(q(tq)), s(ts))
     nb = node_memsize(Q(bq), p)
     nper = tree_n_perm(q(q(tq)), s(ts))
     nsym = mol_block_nsym(q(bq))
     nw = MAX(Hungarian_worksize(n - p, n - p), sdmin_worksize())
!
     if (n == p) then
       call expand_terminal(p, n, nb, nw, nsym, q(q(cq)), q(bq), s(q(ps)), X, X(np), X(nn), W)
     else
       call expand(p, n, nb, nw, nper, nsym, q(q(cq)), q(bq), s(q(ps)), X, X(np), X(nn), W(1), W(2))
     end if
!
     call tree_select_top_node(q(q(tq)), s(ts), UB, X(q(tx)))
!
     if(bb_block_queue_is_empty(q(q(tq)), s(ts))) return
!
     call cswap(p, tree_current_iper(q(q(tq)), s(ts)), s(q(ps)))
!
  end subroutine bb_block_expand
!
  pure subroutine expand(p, n, nb, nw, nper, nsym, cq, b, prm, C, NP, NN, W1, W2)
    integer(IK), intent(in)     :: p, n, nb, nw, nper, nsym
    integer(IK), intent(in)     :: cq(*), b(*), prm(*)
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
      call c_matrix_add(cq, p, prm(iper + p - 1), isym, C, NN(mmap_G, isym, iper), NN(mmap_C, isym, iper))
    end do
!
    do iper = 1, nper
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
    end do
!
  end subroutine expand
!
  pure subroutine expand_terminal(p, n, nb, nw, nsym, cq, b, prm, C, NP, NN, W)
    integer(IK), intent(in)     :: p, n, nb, nw, nsym
    integer(IK), intent(in)     :: cq(*), b(*), prm(*)
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
      call c_matrix_add(cq, p, prm(n), isym, C, NN(mmap_G, isym), NN(mmap_C, isym))
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
!
    do concurrent(j=2:n)
      do concurrent(i=1:r - 1)
        Y(i, j - 1) = X(i, j)
      end do
      do concurrent(i=r + 1:m)
        Y(i - 1, j - 1) = X(i, j)
      end do
    end do
!
  end subroutine subm
!
  pure subroutine cswap(p, iper, prm)
    integer(IK), intent(in)    :: p, iper
    integer(IK), intent(inout) :: prm(*)
    integer(IK)                :: q, i, t
    ! cyclic swap prm(p:q)
    if (iper < 1) return
    q = p + iper
    t = prm(q)
    do i = q, p + 1, -1
      prm(i) = prm(i - 1)
    end do
    prm(p) = t
  end subroutine cswap
!
  pure subroutine creverse(p, iper, prm)
    integer(IK), intent(in)    :: p, iper
    integer(IK), intent(inout) :: prm(*)
    integer(IK)                :: q, i, t
    ! cyclic reverse prm(p:q)
    if (iper < 1) return
    q = p + iper
    t = prm(p)
    do i = p + 1, q
      prm(i - 1) = prm(i)
    end do
    prm(q) = t
!
  end subroutine creverse
!
!| Leave current node.
  pure function bb_block_queue_is_empty(q, s) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! work integer array
    logical                 :: res
    res = tree_queue_is_empty(q(q(tq)), s)
  end function bb_block_queue_is_empty
!
!| Leave current node.
  pure function bb_block_current_value(q, s, X) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! work integer array
    real(RK), intent(in)    :: X(*)
    !! main memory
    real(RK)                :: res
    integer(IK)             :: t
    t = Q(tx) - 1 + tree_current_pointer(q(q(tq)), s(ts))
    res = X(t)
  end function bb_block_current_value
!
!| Leave current node.
  pure function bb_block_queue_is_bottom(q, s) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! work integer array
    logical                 :: res
    res = tree_queue_is_bottom(q(q(tq)), s(ts))
  end function bb_block_queue_is_bottom
!
!| Leave current node.
  pure subroutine bb_block_leave(q, UB, X, s)
    integer(IK), intent(in)    :: q(*)
    !! integer array
    real(RK), intent(in)       :: UB
    !! upperbound
    real(RK), intent(in)       :: X(*)
    !! workarray
    integer(IK), intent(inout) :: s(*)
    !! work integer array
!
    call creverse(tree_current_level(s), tree_current_iper(q(q(tq)), s(ts)), s(q(ps)))
    call tree_leave(q(q(tq)), s(ts))
    call tree_select_top_node(q(q(tq)), s(ts), UB, X(q(tx)))
!
  end subroutine bb_block_leave
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
! update H and C
!     w(it) = HP + H
!     HP = w(it)
!     call add(DD, CP, C, W(ic))
!     call DCOPY(DD, W(ic), 1, CP, 1)
!
! get squared displacement
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
! util
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
!| destructer
  pure elemental subroutine bb_block_destroy(this)
    type(bb_block), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine bb_block_destroy
!
end module mod_bb_block
