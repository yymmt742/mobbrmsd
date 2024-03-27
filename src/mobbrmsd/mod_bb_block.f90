!| Module for handling C, F and tree.
module mod_bb_block
  use mod_params, only: ND, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_params, only: gemm, dot, copy
  use mod_mol_block
  use mod_c_matrix
  use mod_f_matrix
  use mod_rotation
  use mod_Hungarian
  use mod_tree
  implicit none
  private
  public :: bb_block
  public :: bb_block_nmol
  public :: bb_block_molsize
  public :: bb_block_memsize
  public :: bb_block_worksize
  public :: bb_block_setup
  public :: bb_block_inheritance
  public :: bb_block_expand
  public :: bb_block_leave
  public :: bb_block_queue_is_empty
  public :: bb_block_queue_is_bottom
  public :: bb_block_current_value
  public :: bb_block_lowest_value
  public :: bb_block_log_ncomb
  public :: bb_block_evaluation_count
  public :: bb_block_save_state
!
  integer(IK), parameter :: mmap_L = 1
  integer(IK), parameter :: mmap_G = 2
  integer(IK), parameter :: mmap_C = 3
!
  integer(IK), parameter :: header_size = 7
!
  integer(IK), parameter :: cq = 1
  !! pointer to c_matrix interger array
  integer(IK), parameter :: fq = 2
  !! pointer to f_matrix interger array
  integer(IK), parameter :: tq = 3
  !! pointer to tree interger array
  integer(IK), parameter :: fx = 4
  !! pointer to f_matrix memory
  integer(IK), parameter :: tx = 5
  !! pointer to tree memory
  integer(IK), parameter :: cw = 6
  !! pointer to c_matrix work memory
  integer(IK), parameter :: fw = 7
  !! pointer to f_matrix work memory
!
  integer(IK), parameter :: header_memsize = 1
!
  integer(IK), parameter :: nc = 1
  !! pointer to evaluation count (fixed)
  integer(IK), parameter :: cx = 2
  !! pointer to c_matrix array (fixed)
  integer(IK), parameter :: ts = 1
  !! pointer to tree interger work array (fixed)
  integer(IK), parameter :: bq = header_size + 1
  !! pointer to mol_block interger array (fixed)
!
!| bb_block<br>
!  - Node data [L, G, C]<br>
!  - 1    L      : scalar, lowerbound.<br>
!  - 2    G      : scalar, partial sum of auto variance.<br>
!  - 3    C(D,D) : partial sum of covariance.<br>
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
    t = tree(b%q)
!
    q(cq) = bq + SIZE(b%q)
    q(fq) = q(cq) + SIZE(c%q)
    q(tq) = q(fq) + SIZE(f%q)
!
    q(cw) = cx + c_matrix_memsize(c%q)
    q(fx) = q(cw)
    q(fw) = q(fx) + f_matrix_memsize(f%q)
    q(tx) = q(fw)
!
    allocate (res%q, source=[q, b%q, c%q, f%q, t%q])
    allocate (res%s, source=t%s)
!
  end function bb_block_new
!
!| Returns the memory size array size.
  pure function bb_block_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block header array.
    integer(IK)             :: res
!
    res = header_memsize &
   &    + c_matrix_memsize(q(q(cq))) &
   &    + f_matrix_memsize(q(q(fq))) &
   &    + tree_nnodes(q(q(tq))) * ND
!
  end function bb_block_memsize
!
!| Returns the work memory size array size.
  pure function bb_block_worksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block header array.
    integer(IK)             :: p, nmol, nsym, buf, tmp, swrk, hwrk
    integer(IK)             :: res
!
    nmol = mol_block_nmol(q(bq))
    nsym = mol_block_nsym(q(bq))
    swrk = sdmin_worksize()
!
    buf = bb_block_memsize(q)
!
    res = 0
    do p = 1, nmol
      hwrk = Hungarian_worksize(p, p)
      tmp = MAX(MAX(swrk, p**2 + hwrk), swrk * MAX(1, nsym - 1) + 1)
      res = MAX(res, buf + tmp)
      buf = buf - p * nsym * ND
    end do
    res = MAX(res, buf + f_matrix_worksize(q(q(fq))))
    buf = buf - f_matrix_memsize(q(q(fq)))
    res = MAX(res, buf + c_matrix_worksize(q(q(cq))))
    res = res - bb_block_memsize(q)
!
  end function bb_block_worksize
!
!| Returns the number of molecules.
  pure function bb_block_nmol(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block header array.
    integer(IK)             :: res
    res = mol_block_nmol(q(bq))
  end function bb_block_nmol
!
!| Returns the memory size of molecular block size.
  pure function bb_block_molsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block header array.
    integer(IK)             :: res
    res = mol_block_total_size(q(bq))
  end function bb_block_molsize
!
!| Setup C matrix and F matrix in root node.
  pure subroutine bb_block_setup(q, X, Y, s, W, zfill)
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
    logical, intent(in)        :: zfill
    !! if true, the root node is filled by zero.
!
    W(nc) = ZERO
    call tree_reset(q(q(tq)), s(ts))
    call c_matrix_eval(q(q(cq)), q(bq), X, Y, W(cx), W(q(cw)))
    call f_matrix_eval(q(q(fq)), q(q(cq)), W(cx), W(q(fx)), W(q(fw)))
!
    if(.not.zfill) return
!
    block
      real(RK)    :: ZEROS(ND)
      integer(IK) :: nmol
      nmol = mol_block_nmol(q(bq))
      call zfill(ND, ZEROS, 1)
      call evaluate_nodes(nmol, q(q(cq)), q(q(tq)), s(ts), W(cx), W(q(fx)), ZEROS, W(q(tx)), W(nc))
      call tree_select_top_node(q(q(tq)), s(ts), ND, RHUGE, W(q(tx)))
    end block
!
  end subroutine bb_block_setup
!
!| Expands the latest node of the parent block to the top-level queue of the child block.
  pure subroutine bb_block_inheritance(UB, q, s, X, p, r, Z)
    real(RK), intent(in)       :: UB
    !! upper bound
    integer(IK), intent(in)    :: q(*)
    !! work integer array
    integer(IK), intent(inout) :: s(*)
    !! work integer array
    real(RK), intent(inout)    :: X(*)
    !! main memory
    integer(IK), intent(in)    :: p(*)
    !! parent integer array
    integer(ik), intent(in)    :: r(*)
    !! parent integer work array
    real(RK), intent(in)       :: Z(*)
    !! parent work array
    integer(IK)                :: pp, nmol
!
    nmol = mol_block_nmol(q(bq))
    pp = p(tx) + (tree_current_pointer(p(p(tq)), r(ts)) - 1) * ND
!
    call evaluate_nodes(nmol, q(q(cq)), q(q(tq)), s(ts), X(cx), X(q(fx)), Z(pp), X(q(tx)), X(nc))
    call tree_reset(q(q(tq)), s(ts))
    call tree_select_top_node(q(q(tq)), s(ts), ND, UB, X(q(tx)))
!
  end subroutine bb_block_inheritance
!
!| Expand top node in queue.
  pure subroutine bb_block_expand(UB, q, s, X)
    real(RK), intent(in)       :: UB
    !! upper bound
    integer(IK), intent(in)    :: q(*)
    !! work integer array
    integer(IK), intent(inout) :: s(*)
    !! work integer array
    real(RK), intent(inout)    :: X(*)
    !! main memory
    integer(IK)                :: nmol, pp
!
     nmol = mol_block_nmol(q(bq))
     block
       do
         if (bb_block_queue_is_empty(q, s) .or. bb_block_queue_is_bottom(q, s)) return
         pp = q(tx) + (tree_current_pointer(q(q(tq)), s(ts)) - 1) * ND
         call tree_expand(q(q(tq)), s(ts))
         call evaluate_nodes(nmol, q(q(cq)), q(q(tq)), s(ts), X(cx), X(q(fx)), X(pp), X(q(tx)), X(nc))
         call tree_select_top_node(q(q(tq)), s(ts), ND, UB, X(q(tx)))
       end do
     end block
!
  end subroutine bb_block_expand
!
!| leave current node.
  pure subroutine bb_block_leave(UB, q, s, X)
  !pure subroutine bb_block_leave(UB, q, s, X)
    real(RK), intent(in)       :: UB
    !! upper bound
    integer(IK), intent(in)    :: q(*)
    !! work integer array
    integer(IK), intent(inout) :: s(*)
    !! work integer array
    real(RK), intent(in)       :: X(*)
    !! main memory
    integer(IK)                :: l
!
     do
       call tree_select_top_node(q(q(tq)), s(ts), ND, UB, X(q(tx)))
       l = tree_current_level(s(ts))
       if (l == 1 .or. .not. bb_block_queue_is_empty(q, s)) return
       call tree_leave(q(q(tq)), s(ts))
     enddo
!
  end subroutine bb_block_leave
!
  pure subroutine evaluate_nodes(nmol, qcov, q, s, CX, FX, Z, X, neval)
    integer(IK), intent(in) :: nmol, qcov(*), q(*), s(*)
    real(RK), intent(in)    :: CX(*), FX(*), Z(ND)
    real(RK), intent(inout) :: X(ND, *), neval
    integer(IK)             :: l, px, pw, m, nw, nper, nsym, iper, perm(nmol)
!
    l = tree_current_level(s)
    m = nmol - l
    nw = sdmin_worksize()
    px = tree_queue_pointer(q, s)
!
    nsym = tree_n_sym(q)
    nper = tree_n_perm(q, s)
    perm = tree_current_permutation(q, s)
!
    pw = px + nsym
!
    do iper = 1, nper
      call evaluate_queue(iper, l, m, nw, nsym, nper, nmol, qcov, perm, &
     &                    CX, FX, Z, X(1, px), X(1, pw), X(2, pw))
      pw = pw + nsym
    end do
!
    neval = neval + real(nsym * nper, RK)
!
  end subroutine evaluate_nodes
!
  pure subroutine evaluate_queue(iper, l, m, nw, nsym, nper, nmol, qcov, perm, CX, FX, Z, X, W1, W2)
    integer(IK), intent(in) :: iper, l, m, nw, nsym, nper, nmol
    integer(IK), intent(in) :: qcov(*), perm(*)
    real(RK), intent(in)    :: CX(*), FX(*), Z(ND)
    real(RK), intent(inout) :: X(ND, nsym, nper), W1(*), W2(nw, *)
    integer(IK)             :: isym, iprm
!
    iprm = perm(l + iper - 1)
!
    call copy(ND, Z, 1, X(1, 1, iper), 1)
    call c_matrix_add(qcov, l, iprm, 1, CX, X(mmap_G, 1, iper), X(mmap_C, 1, iper))
    call subm(l, nmol, iper, perm, FX, W1)
    call Hungarian(m, m, W1(1), W1(m * m + 1))
!
    do concurrent(isym=2:nsym)
      call copy(ND, Z, 1, X(1, isym, iper), 1)
      call c_matrix_add(qcov, l, iprm, isym, CX, X(mmap_G, isym, iper), X(mmap_C, isym, iper))
      call estimate_sdmin(X(mmap_G, isym, iper), X(mmap_C, isym, iper), W2(1, isym - 1))
      X(mmap_L, isym, iper) = W1(1) + W2(1, isym - 1)
    end do
!
    call estimate_sdmin(X(mmap_G, 1, iper), X(mmap_C, 1, iper), W2)
    X(mmap_L, 1, iper) = W1(1) + W2(1, 1)
!
  end subroutine evaluate_queue
!
  pure subroutine subm(l, n, r, perm, X, Y)
    integer(IK), intent(in) :: l, n, r, perm(*)
    real(RK), intent(in)    :: X(n, n)
    real(RK), intent(inout) :: Y(n - l, n - l)
    integer(IK)             :: i, j
!
    do concurrent(j=l + 1:n)
      do concurrent(i=l:l + r - 2)
        Y(i - l + 1, j - l) = X(perm(i), j)
      end do
      do concurrent(i=l + r:n)
        Y(i - l, j - l) = X(perm(i), j)
      end do
    end do
!
  end subroutine subm
!
!| Returns true when queue is empty
  pure function bb_block_queue_is_empty(q, s) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! work integer array
    logical                 :: res
    res = tree_queue_is_empty(q(q(tq)), s(ts))
  end function bb_block_queue_is_empty
!
!| Returns true when queue is bottom
  pure function bb_block_queue_is_bottom(q, s) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! work integer array
    logical                 :: res
    res = tree_queue_is_bottom(q(q(tq)), s(ts))
  end function bb_block_queue_is_bottom
!
!| Returns current_value
  pure function bb_block_current_value(q, s, X) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! work integer array
    real(RK), intent(in)    :: X(*)
    !! main memory
    real(RK)                :: res
    integer(IK)             :: t
    t = q(tx) + ND * (tree_current_pointer(q(q(tq)), s(ts)) - 1)
    res = X(t)
  end function bb_block_current_value
!
!| Returns the minimum value of the surviving nodes, excluding the current value.
  pure function bb_block_lowest_value(q, s, X) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! state vector
    real(RK), intent(in)    :: X(*)
    !! main memory
    real(RK)                :: res
    res = tree_lowest_value(q(q(tq)), s(ts), ND, X(q(tx)))
  end function bb_block_lowest_value
!
!| Returns the minimum value of the surviving nodes, excluding the current value.
  pure function bb_block_evaluation_count(X) result(res)
    real(RK), intent(in)    :: X(*)
    !! main memory
    real(RK)                :: res
    res = X(nc)
  end function bb_block_evaluation_count
!
!| Returns the minimum value of the surviving nodes, excluding the current value.
  pure function bb_block_log_ncomb(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    real(RK)                :: res
    res = tree_log_ncomb(q(q(tq)))
  end function bb_block_log_ncomb
!
!| Save current state.
  pure subroutine bb_block_save_state(q, s, z)
    integer(IK), intent(in)    :: q(*)
    !! integer array
    integer(IK), intent(in)    :: s(*)
    !! state vector
    integer(IK), intent(inout) :: z(*)
    !! memory
    integer(IK)                :: nmol
    nmol = mol_block_nmol(q(bq))
    z(:nmol) = tree_current_sequence(q(q(tq)), s(ts))
  end subroutine bb_block_save_state
!
! util
!
!| destructer
  pure elemental subroutine bb_block_destroy(this)
    type(bb_block), intent(inout) :: this
    if (ALLOCATED(this%q)) deallocate (this%q)
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine bb_block_destroy
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
end module mod_bb_block

