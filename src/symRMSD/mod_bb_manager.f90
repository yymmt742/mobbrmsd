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
  use mod_lowerbound
  use mod_tree
  implicit none
  private
  public :: bb_manager
  public :: bb_manager_tuple
  public :: bb_manager_memsize
  public :: bb_manager_setup
  public :: bb_manager_setup_root
  public :: bb_manager_expand
  public :: bb_manager_worksize
!
!| bb_manager<br>
!    Node data [L, G, C, F, W]<br>
!    0    L      : scalar, lowerbound.<br>
!    1    G      : scalar, partial sum of auto variance.<br>
!    2    F(n,n) : Free rotate score matrix.<br>
!    2+nn C(D,D) : partial sum of covariance.<br>
!         W(*)   : work space.<br>
  type bb_manager
    sequence
    private
    integer(IK), public :: x
    !! pointer offset of rwork.
    integer(IK), public :: q
    !! pointer offset of iwork.
    integer(IK)         :: nm
    !! memory size.
    integer(IK)         :: nw
    !! work size.
    type(mol_block)     :: b
    !! mol_block
    type(c_matrix)      :: c
    !! c_matrix
    type(f_matrix)      :: f
    !! f_matrix
    type(tree)          :: t
    !! tree
  end type bb_manager
!
!| A set of bb_manager and work arrays. <br>
!  This is mainly used for passing during initialization.
  type bb_manager_tuple
    type(bb_manager)         :: bb
    !! bb_manager
    integer(IK), allocatable :: q(:)
    !! work integer array
    real(RK), allocatable    :: x(:)
    !! main memory
  contains
    final           :: bb_manager_tuple_destroy
  end type bb_manager_tuple
!
  interface bb_manager_tuple
    module procedure bb_manager_tuple_new
  end interface bb_manager_tuple
!
contains
!
!| Constructer
  pure function bb_manager_tuple_new(m, n, sym) result(res)
    integer(IK), intent(in)           :: m
    !! number of molecules.
    integer(IK), intent(in)           :: n
    !! number of atoms per molecule.
    integer(IK), intent(in), optional :: sym(:,:)
    !! symmetric codomains, [[a1,a2,...,am], [b1,b2,...,bm], ...].
    type(mol_block_tuple)             :: b
    type(c_matrix_tuple)              :: c
    type(f_matrix_tuple)              :: f
    type(tree_tuple)                  :: t
    type(bb_manager_tuple)            :: res
!
    b = mol_block_tuple(m, n, sym)
    c = c_matrix_tuple(b%b)
    f = f_matrix_tuple(b%b)
    t = tree_tuple(b%b, node_memsize)
!
    res%bb%x = 1
    res%bb%q = 1
!
    res%bb%b = b%b
    res%bb%c = c%c
    res%bb%f = f%f
    res%bb%t = t%t
!
    res%bb%nm = SIZE(c%x) + SIZE(f%x) + SIZE(t%x)
    res%bb%nw = SIZE(f%x)
    res%bb%nw = MAX(res%bb%nw, SIZE(c%x) - res%bb%nw)
    res%bb%nw = MAX(res%bb%nw - SIZE(t%x), 0)
!
    res%bb%t%q = SIZE(b%w) + 1
!
    res%bb%c%p = 1
    res%bb%c%w = res%bb%c%p + c_matrix_memsize(res%bb%c)
    res%bb%f%p = res%bb%c%w
    res%bb%f%w = res%bb%f%p + f_matrix_memsize(res%bb%f)
    res%bb%t%x = res%bb%f%w
!
    allocate (res%q, source=[b%w, t%q])
    allocate (res%x(res%bb%nm + res%bb%nw))
!
  end function bb_manager_tuple_new
!
  pure function node_memsize(b, p) result(res)
    type(mol_block), intent(in) :: b
    integer(IK), intent(in)     :: p
    integer(IK)                 :: res, n
    n = mol_block_nmol(b) - p
    res = 1 + n * n + 1 + DD
  end function node_memsize
!
!| Inquire worksize of f_matrix.
  pure elemental function bb_manager_memsize(this) result(res)
    type(bb_manager), intent(in) :: this
    !! bb_manager.
    integer(IK)                  :: res
    res = this%nm
  end function bb_manager_memsize
!
!| Inquire worksize of bb_manager.
  pure elemental function bb_manager_worksize(this) result(res)
    type(bb_manager), intent(in) :: this
    !! bb_manager.
    integer(IK)                  :: res
    res = this%nw
  end function bb_manager_worksize
!
  pure function root_pointer(this, Q) result(res)
    type(bb_manager), intent(in) :: this
    integer(IK), intent(in)      :: Q(*)
    integer(IK)                  :: res
    res = this%x + tree_root_pointer(this%t, Q(this%q))
  end function root_pointer
!
!| Setup C matrix and F matrix in root node.
  pure subroutine bb_manager_setup(this, Q, X, Y, W)
    type(bb_manager), intent(in) :: this
    !! bb_manager
    integer(IK), intent(in)      :: Q(*)
    !! work integer array
    real(RK), intent(in)         :: X(*)
    !! reference coordinate
    real(RK), intent(in)         :: Y(*)
    !! target coordinate
    real(RK), intent(inout)      :: W(*)
    !! work integer array
    integer(IK)                  :: ix
!
    ix = mol_block_pointer(this%b)
    call c_matrix_eval(this%c, this%b, Q(this%q), X(ix), Y(ix), W(this%x), W(this%x))
    call f_matrix_eval(this%f, this%b, this%c, W(this%x), W(this%x), W(this%x))
!
    ix = root_pointer(this, Q)
    W(ix) = ZERO
    ix = ix + 2 + DD
    call copy(mol_block_nmol(this%b)**2, W(this%x + this%f%p - 1), 1, W(ix), 1)
!
  end subroutine bb_manager_setup
!
!| Setup C matrix and F matrix in root node.
  pure subroutine bb_manager_setup_root(this, Q, W, G0, C0)
    type(bb_manager), intent(in)   :: this
    !! bb_manager
    integer(IK), intent(in)        :: Q(*)
    !! work integer array
    real(RK), intent(inout)        :: W(*)
    !! work array
    real(RK), intent(in), optional :: G0
    !! partial sum of auto variances.
    real(RK), intent(in), optional :: C0(*)
    !! main memory
    integer(IK)                    :: t
!
    t = root_pointer(this, Q) + 1
!
    if (PRESENT(G0)) then; W(t) = G0
    else; W(t) = ZERO
    end if
!
    if (PRESENT(C0)) then; call copy(DD, C0, 1, W(t + 1), 1)
    else; call zfill(DD, W(t + 1), 1)
    end if
!
  end subroutine bb_manager_setup_root
!
!| Expand top node in queue.
  pure subroutine bb_manager_expand(this, Q, W)
    type(bb_manager), intent(in) :: this
    !! bb_manager
    integer(IK), intent(inout)   :: Q(*)
    !! work integer array
    real(RK), intent(inout)      :: W(*)
    !! work real array
    integer(IK)                  :: nper, nsym
    integer(IK)                  :: p, n, np, nn, nb
!
     np = this%x + tree_current_pointer(this%t, Q(this%q))
     call tree_expand(this%t, Q(this%q))
     nn = this%x + tree_queue_pointer(this%t, Q(this%q))
!
     p = tree_current_level(this%t, Q(this%q))
     n = mol_block_nmol(this%b)
     nb = node_memsize(this%b, p)
     nper = tree_n_perm(this%t, Q(this%q))
     nsym = mol_block_nsym(this%b)
!
     call expand(this%c, this%b, p, n, nb, nper, nsym, Q(this%q), W(this%x), W(np), W(nn))
!
  end subroutine bb_manager_expand
!
  pure subroutine expand(cm, b, p, n, nb, nper, nsym, Q, C, NP, NN)
    type(c_matrix), intent(in)  :: cm
    type(mol_block), intent(in) :: b
    integer(IK), intent(in)     :: p, n, nb, nper, nsym
    integer(IK), intent(in)     :: Q(*)
    real(RK), intent(in)        :: C(*)
    real(RK), intent(in)        :: NP(*)
    real(RK), intent(inout)     :: NN(nb, nsym, nper)
    integer(IK), parameter      :: mmap_L = 1
    integer(IK), parameter      :: mmap_G = 2
    integer(IK), parameter      :: mmap_C = 3
    integer(IK)                 :: mmap_F
    integer(IK)                 :: mp, mn, mm, nw1, nw2
    integer(IK)                 :: iper, isym
!
    mn = n - p
    mm = mn**2
    mmap_F = 3 + DD
    mp = mn + 1
    nw1 = MAX(Hungarian_worksize(mn, mn), sdmin_worksize())
    nw2 = sdmin_worksize()
!
    do concurrent(iper=1:nper)
      block
        real(RK) :: w1(nw1)
!
        call subm(mp, mp, iper, NP(mmap_F), NN(mmap_F, 1, iper))
        call Hungarian(mn, mn, NN(mmap_F, 1, iper), w1)
        NN(mmap_L, 1, iper) = w1(1)
!
        call copy(DD + 1, NP(mmap_G), 1, NN(mmap_G, 1, iper), 1)
        call c_matrix_add(cm, b, iper, p, 1, C, NN(mmap_G, 1, iper), NN(mmap_C, 1, iper))
        call estimate_sdmin(NN(mmap_G, 1, iper), NN(mmap_C, 1, iper), w1)
!
        do concurrent(isym=2:nsym)
!
          block
            real(RK) :: w2(nw2)
!
            call copy(DD + 1, NP(mmap_G), 1, NN(mmap_G, isym, iper), 1)
            call copy(mm, NN(mmap_F, 1, iper), 1, NN(mmap_F, isym, iper), 1)
!
            call c_matrix_add(cm, b, iper, p, isym, C, NN(mmap_G, isym, iper), NN(mmap_C, isym, iper))
            call estimate_sdmin(NN(mmap_G, isym, iper), NN(mmap_C, isym, iper), w2)
            NN(mmap_L, isym, iper) = w2(1) + NN(mmap_L, 1, iper)
!
          end block
        end do
!
        NN(mmap_L, 1, iper) = w1(1) + NN(mmap_L, 1, iper)
!
      end block
    end do
!
!
  end subroutine expand
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
!| lowerbound function.
!  L(G, C, D) = SUM_{i=1,...,p} (G - 2tr[CR]) + min_{nu} SUM_{i=p+1,...,N} D_{i nu(p)}
! pure subroutine lowerbound(p, b, G, C, D, W)
!   integer(IK), intent(in)    :: p
!   !! level
!   type(mol_block),intent(in) :: b
!   !! mol_block
!   real(RK), intent(in)       :: G
!   !! partial auto variance, G
!   real(RK), intent(in)       :: C(*)
!   !! partial covariance matrix, C(d, d)
!   real(RK), intent(in)       :: D(*)
!   !! residual matrix, D(n1, n2), here n1 = MAX(nx, ny) - p and n2 = MIN(nx, ny) - p.
!   real(RK), intent(inout)    :: W(*)
!   !! workarray, must be SIZE(W) > lowerbound_worksize(p, b).
!   integer(IK)                :: n
!
!   W(1) = ZERO
!   n = mol_block_nmol(b) - p
!
!   if (p < 0 .or. n < 0) return
!   if (0 < n) call Hungarian(n, n, D, W)
!   if (0 < p) then
!     call estimate_sdmin(G, C, W(2))
!     W(1) = W(1) + W(2)
!   end if
!
! end subroutine lowerbound
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
  pure elemental subroutine bb_manager_tuple_destroy(this)
    type(bb_manager_tuple), intent(inout) :: this
    if (ALLOCATED(this%x)) deallocate (this%x)
    if (ALLOCATED(this%q)) deallocate (this%q)
  end subroutine bb_manager_tuple_destroy
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

