!| Module for handling \( \{\mathbf{C}\}_{IJs}, \mathbf{F} \) and the factrial tree. <br>
!  This module provides procedures for performing branch-and-bound algorithms (BB)
!  on assemblies consisting of homologous molecules. <br>
!  Since the BB procedure is implemented for multi-component systems,
!  this module does not provide the BB itself. <br>
!  @note
!    Node data block is defined by \( [l_p, G_p, \mathbf{C}_p]\), here <br>
!    \( l_p \)          : scalar, lowerbound.<br>
!    \[ l_p = - \max_{\mathbf{R}} \text{Tr} \left[ \mathbf{C}_p\mathbf{R}\right] - \max_{\nu'}\sum_{I=p+1}^M F_{I\nu'(I)} \]
!    \( G_p \)          : scalar, partial sum of auto variance.<br>
!    \[ G_p = \sum_{I=1}^p \mathbf{G}_{I\nu(I)} \]
!    \( \mathbf{C}_p \) : partial sum of covariance. <br>
!    \[ \mathbf{C}_p = \sum_{I=1}^p \mathbf{C}_{I\nu(I)\sigma(I)} \]
!  @endnote
module mod_bb_block
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use blas_lapack_interface, only: D, ND
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
  public :: bb_block_natm
  public :: bb_block_molsize
  public :: bb_block_memsize
  public :: bb_block_worksize
  public :: bb_block_setup
  public :: bb_block_inheritance
  public :: bb_block_expand
  public :: bb_block_closure
  public :: bb_block_tree_is_empty
  public :: bb_block_tree_is_bottom
  public :: bb_block_current_value
  public :: bb_block_lowest_value
  public :: bb_block_log_ncomb
  public :: bb_block_evaluation_count
  public :: bb_block_save_state
  public :: bb_block_swap_y
  public :: bb_block_covmat_add
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
  integer(IK), parameter :: wcov = 2
  !! pointer to c_matrix array (fixed)
  integer(IK), parameter :: stree = 1
  !! pointer to tree interger work array (fixed)
  integer(IK), parameter :: qblck = header_size + 1
!
!| bb_block<br>
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
  pure function bb_block_new(n, M, sym) result(res)
    integer(IK), intent(in) :: n
    !! number of molecules.
    integer(IK), intent(in) :: M
    !! number of atoms per molecule.
    integer(IK), intent(in), optional :: sym(:, :)
    !! symmetric codomains, [[a1,a2,...,am], [b1,b2,...,bm], ...].
    type(mol_block)         :: b
    type(c_matrix)          :: c
    type(f_matrix)          :: f
    type(tree)              :: t
    type(bb_block)          :: res
    integer(IK)             :: q(header_size)
!
    b = mol_block(n, M, sym)
    c = c_matrix(b%q)
    f = f_matrix(b%q)
    t = tree(mol_block_nmol(b%q), mol_block_nsym(b%q))
!
    q(cq) = qblck + SIZE(b%q)
    q(fq) = q(cq) + SIZE(c%q)
    q(tq) = q(fq) + SIZE(f%q)
!
    q(cw) = wcov + c_matrix_memsize(c%q)
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
    nmol = mol_block_nmol(q(qblck))
    nsym = mol_block_nsym(q(qblck))
    swrk = sdmin_worksize()
!
    buf = bb_block_memsize(q)
!
    res = 0
    do p = 1, nmol
      hwrk = Hungarian_worksize(p, p)
      tmp = MAX(MAX(swrk, p**2 + hwrk), swrk * MAX(1, nsym) + 1)
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
    res = mol_block_nmol(q(qblck))
  end function bb_block_nmol
!
!| Returns the number of molecules.
  pure function bb_block_natm(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block header array.
    integer(IK)             :: res
    res = mol_block_natm(q(qblck))
  end function bb_block_natm
!
!| Returns the memory size of molecular block size.
  pure function bb_block_molsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! bb_block header array.
    integer(IK)             :: res
    res = mol_block_total_size(q(qblck))
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
    call tree_reset(q(q(tq)), s(stree))
    call c_matrix_eval(q(q(cq)), q(qblck), X, Y, W(wcov), W(q(cw)))
    call f_matrix_eval(q(q(fq)), q(q(cq)), W(wcov), W(q(fx)), W(q(fw)))
!
    if (.not. zfill) return
!
    block
      integer(IK) :: nmol
      real(RK)    :: ZEROS(ND)
      ZEROS = ZERO
      nmol = mol_block_nmol(q(qblck))
      call evaluate_nodes(nmol, q(q(cq)), q(q(tq)), s(stree), W(wcov), W(q(fx)), ZEROS, W(q(tx)), W(nc))
!     call tree_select_top_node(q(q(tq)), s(stree), ND, RHUGE, W(q(tx)))
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
    nmol = mol_block_nmol(q(qblck))
    pp = p(tx) + (tree_current_pointer(p(p(tq)), r(stree)) - 1) * ND
!
    call evaluate_nodes(nmol, q(q(cq)), q(q(tq)), s(stree), X(wcov), X(q(fx)), Z(pp), X(q(tx)), X(nc))
    call tree_reset(q(q(tq)), s(stree))
!   call tree_select_top_node(q(q(tq)), s(stree), ND, UB, X(q(tx)))
!
  end subroutine bb_block_inheritance
!
!| Expand top node in queue.
  pure subroutine bb_block_expand(UB, q, s, W)
    real(RK), intent(in)       :: UB
    !! upper bound
    integer(IK), intent(in)    :: q(*)
    !! header
    integer(IK), intent(inout) :: s(*)
    !! state
    real(RK), intent(inout)    :: W(*)
    !! main memory
    integer(IK)                :: nmol
    integer(IK)                :: pp ! previous current node
    associate ( &
   &  qtree => q(tq),&
   &  wtree => q(tx),&
   &  qcov => q(cq), &
   &  wfrx => q(fx)&
   &  )
      nmol = mol_block_nmol(q(qblck))
      block
        do
          call tree_select_top_node(q(q(tq)), s(stree), ND, UB, W(q(tx)))
          if (tree_queue_is_bottom(q(qtree), s(stree)) &
       & .or. tree_queue_is_empty(q(qtree), s(stree))) exit
          pp = q(tx) + (tree_current_pointer(q(qtree), s(stree)) - 1) * ND
          call tree_expand(q(qtree), s(stree))
          call evaluate_nodes(nmol, q(qcov), q(qtree), s(stree), W(wcov), W(wfrx), W(pp), W(wtree), W(nc))
        end do
      end block
    end associate
  end subroutine bb_block_expand
!
!| closure current node.
  pure subroutine bb_block_closure(UB, q, s, X)
    real(RK), intent(in)       :: UB
    !! upper bound
    integer(IK), intent(in)    :: q(*)
    !! work integer array
    integer(IK), intent(inout) :: s(*)
    !! work integer array
    real(RK), intent(in)       :: X(*)
    !! main memory
    associate (qtree => q(tq), wtree => q(tx))
      do
        if (tree_queue_is_left(q(qtree), s(stree)) &
     & .or. tree_queue_is_root(q(qtree), s(stree))) exit
        call tree_leave(q(qtree), s(stree))
      end do
    end associate
  end subroutine bb_block_closure
!
!| Evaluate nodes in tree_current_level.
  pure subroutine evaluate_nodes(nmol, qcov, q, s, C, F, Z, W, neval)
    integer(IK), intent(in) :: nmol, qcov(*), q(*), s(*)
    real(RK), intent(in)    :: C(*), F(*), Z(ND)
    real(RK), intent(inout) :: W(ND, *), neval
    integer(IK)             :: l, px, pw, m, nw, nh, nper, nsym, iper, perm(nmol)
    l = tree_current_level(s)
    m = nmol - l ! residual dimension
    nw = sdmin_worksize()
    nh = Hungarian_worksize(m, m)
    px = tree_queue_pointer(q, s)
!
    nsym = tree_n_sym(q)
    nper = tree_n_perm(q, s)
    perm = tree_current_permutation(q, s)
    pw = px + nsym
!
    do iper = 1, nper
      call evaluate_queue(iper, perm(l + iper - 1), l, &
     &                    m, nw, nh, nsym, nper, nmol, qcov, perm, &
     &                    C, F, Z, W(1, px), W(1, pw), W(2, pw))
      pw = pw + nsym
    end do
!
    neval = neval + real(nsym * nper, RK)
  end subroutine evaluate_nodes
!
!| Evaluate nodes in tree_current_level and iper.
  pure subroutine evaluate_queue(iper, iabp, l, m, nw, nh, nsym, nper, nmol, qcov, perm, C, F, Z, X, W1, W2)
    integer(IK), intent(in) :: iper, iabp, l, m, nw, nh, nsym, nper, nmol
    integer(IK), intent(in) :: qcov(*), perm(*)
    real(RK), intent(in)    :: C(*), F(*), Z(ND)
    real(RK), intent(inout) :: X(ND, nsym, nper), W1(*), W2(nw, *)
    integer(IK)             :: isym
    if (m == 0) then
      W1(1) = ZERO
    else
      call subm(l, nmol, iper, perm, F, W1(nh + 1))
      call Hungarian(m, m, W1(nh + 1), W1(1))
    end if
    do concurrent(isym=1:nsym)
      call copy(ND, Z, X(1, isym, iper))
      call c_matrix_add(qcov, l, iabp, isym, C, X(mmap_G, isym, iper), X(mmap_C, isym, iper))
      call estimate_rcmax(X(mmap_G, isym, iper), X(mmap_C, isym, iper), W2(1, isym))
      X(mmap_L, isym, iper) = W1(1) - W2(1, isym)
    end do
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
!| Returns true when tree is empty
  pure function bb_block_tree_is_empty(q, s) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! work integer array
    logical                 :: res
    res = tree_queue_is_empty(q(q(tq)), s(stree)) &
   &.and. tree_queue_is_root(q(q(tq)), s(stree))
  end function bb_block_tree_is_empty
!
!| Returns true when tree is bottom
  pure function bb_block_tree_is_bottom(q, s) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! work integer array
    logical                 :: res
    res = tree_queue_is_left(q(q(tq)), s(stree)) &
   &.and. tree_queue_is_bottom(q(q(tq)), s(stree))
  end function bb_block_tree_is_bottom
!
!| Returns current L value.
  pure function bb_block_current_value(q, s, W) result(res)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: s(*)
    !! work integer array
    real(RK), intent(in)    :: W(*)
    !! main memory
    real(RK)                :: res
    integer(IK)             :: t
    t = q(tx) + mmap_L - 1 + ND * (tree_current_pointer(q(q(tq)), s(stree)) - 1)
    res = W(t)
  end function bb_block_current_value
!
!| Returns the minimum value of the surviving nodes, excluding the current value.
  pure function bb_block_lowest_value(q, s, W) result(res)
    integer(IK), intent(in) :: q(*)
    !! header
    integer(IK), intent(in) :: s(*)
    !! state
    real(RK), intent(in)    :: W(*)
    !! main memory
    real(RK)                :: res
    associate (qtree => q(tq), wtree => q(tx))
      res = tree_lowest_value(q(qtree), s(stree), ND, W(wtree))
    end associate
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
!
    nmol = mol_block_nmol(q(qblck))
    z(:nmol) = tree_current_sequence(q(q(tq)), s(stree))
!
  end subroutine bb_block_save_state
!
!| swap Y by saved state z.
  pure subroutine bb_block_swap_y(q, z, Y)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: z(*)
    !! saved state (not state vector)
    real(RK), intent(inout) :: Y(*)
    !! target coordinate
    integer(IK)             :: nmol, napm
!
    nmol = mol_block_nmol(q(qblck))
    napm = mol_block_napm(q(qblck))
    block
      integer(IK) :: iper(nmol), imap(nmol)
      iper = tree_sequence_to_permutation(q(q(tq)), z)
      imap = tree_sequence_to_mapping(q(q(tq)), z)
      call swap_y(nmol, napm, iper, imap, q(qblck), Y)
    end block
!
  contains
    pure subroutine swap_y(nmol, napm, iper, imap, b, Y)
      integer(IK), intent(in) :: nmol, napm, iper(nmol), imap(nmol), b(*)
      real(RK), intent(inout) :: Y(D, napm, nmol)
      real(RK)                :: T(D, napm, nmol)
      integer(IK)             :: i, dm, dmn
      dm = D * napm
      dmn = dm * nmol
      do concurrent(i=1:nmol)
        call copy(dm, Y(:, :, i), T(:, :, iper(i)))
        call mol_block_swap(b, imap(i), T(1, 1, iper(i)))
      end do
      call copy(dmn, T, Y)
    end subroutine swap_y
  end subroutine bb_block_swap_y
!
!| Sum covariance matrix by saved state z.
  pure subroutine bb_block_covmat_add(q, z, W, G, C)
    integer(IK), intent(in) :: q(*)
    !! integer array
    integer(IK), intent(in) :: z(*)
    !! saved state (not state vector)
    real(RK), intent(in)    :: W(*)
    !! main memory
    real(RK), intent(inout) :: G
    !! autovariance
    real(RK), intent(inout) :: C(*)
    !! covariance matrix
    integer(IK)             :: nmol
    nmol = mol_block_nmol(q(qblck))
    block
      integer(IK) :: i, iper(nmol), imap(nmol)
      iper = tree_sequence_to_permutation(q(q(tq)), z)
      imap = tree_sequence_to_mapping(q(q(tq)), z) + 1
      do i = 1, nmol
        call c_matrix_add(q(q(cq)), i, iper(i), imap(i), W(wcov), G, C)
      end do
    end block
  end subroutine bb_block_covmat_add
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
  pure subroutine copy(N, X, Y)
    integer(IK), intent(in) :: N
    real(RK), intent(in)    :: X(*)
    real(RK), intent(inout) :: Y(*)
    integer(IK)             :: i
    do concurrent(i=1:N)
      Y(i) = X(i)
    end do
  end subroutine copy
!
end module mod_bb_block

