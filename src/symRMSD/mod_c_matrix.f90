!
!| Module for manage C matrix.<br>
!  {C_IJs} :: Covariance matrices.<br>
!  - C_IJs(d,d) = Y_J @ Q_s @ X_I^T<br>
!  - - X_I :: I-th molecule in X.<br>
!  - - Y_J :: J-th molecule in Y.<br>
!  - - Q_s :: Permutation matrix on m.
module mod_c_matrix
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE, &
    &                   gemm, dot, copy, axpy
  use mod_mol_block
  implicit none
  private
  public :: c_matrix
  public :: c_matrix_tuple
  public :: c_matrix_memsize
  public :: c_matrix_worksize
  public :: c_matrix_blocksize
  public :: c_matrix_eval
  public :: c_matrix_add
!
!| c_matrix <br>
!  - C_IJs is the d x d matrix and {C_IJs} is the third-order tensor of nx x ny x s. <br>
!  - To quickly find the rotation matrix, C is stored with G_IJ given by <br>
!  --- G_IJ = Tr[X_I @ X_I^T] + Tr[Y_J @ Y_J^T]. <br>
!  - G does not change with respect to s. <br>
!  - C(:,I,J) = [G_IJ, C_IJ1, C_IJ2, ..., C_IJS] with C_IJs(D,D) := Y_J @ Q_s @ X_I^T. <br>
  type c_matrix
    private
    sequence
    integer(IK), public :: p
    !! pointer to main memory.
    integer(IK), public :: w
    !! pointer to work array.
    integer(IK)         :: nl
    !! number of row. nl = n.
    integer(IK)         :: cb
    !! number of elements in a sell. cb = DD * res%b%s + 1.
    integer(IK)         :: cl
    !! number of elements in a line. cl = cb * n.
    integer(IK)         :: nw
    !! number of work array. nw = MAX(dmn + dm, n + n)
  end type c_matrix
!
  interface c_matrix
    module procedure c_matrix_new
  end interface c_matrix
!
!| A set of c_matrix and work arrays. <br>
!  This is mainly used for passing during initialization.
  type c_matrix_tuple
    type(c_matrix)        :: c
    !! header
    real(RK), allocatable :: x(:)
    !! main memory.
    real(RK), allocatable :: w(:)
    !! work memory.
  end type c_matrix_tuple
!
  interface c_matrix_tuple
    module procedure c_matrix_tuple_new
  end interface c_matrix_tuple
!
contains
!
!| Constructer
  pure elemental function c_matrix_tuple_new(b) result(res)
    type(mol_block), intent(in) :: b
    !! mol_block, must be initialized.
    type(c_matrix_tuple)        :: res
!
    res%c = c_matrix_new(b)
    allocate (res%x(c_matrix_memsize(res%c)))
    allocate (res%w(c_matrix_worksize(res%c)))
!
  end function c_matrix_tuple_new
!
!| Constructer
  pure elemental function c_matrix_new(b) result(res)
    type(mol_block), intent(in) :: b
    !! mol_block, must be initialized.
    type(c_matrix)              :: res
!
    res%p = 1
    res%w = 1
    res%cb = 1 + DD * mol_block_nsym(b)
    res%nl = mol_block_nmol(b)
    res%cl = res%cb * res%nl
    res%nw = MAX(mol_block_each_size(b) + mol_block_total_size(b), &
           &     res%nl + res%nl)
!
  end function c_matrix_new
!
!| Inquire blocksize of c_matrix.
  pure elemental function c_matrix_blocksize(this) result(res)
    type(c_matrix), intent(in) :: this
    !! c_matrix
    integer(IK)                :: res
    res = this%cb
  end function c_matrix_blocksize
!
!| Inquire memsize of c_matrix.
  pure elemental function c_matrix_memsize(this) result(res)
    !| this :: c_matrix
    type(c_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%cb * this%nl * this%nl
  end function c_matrix_memsize
!
!| Inquire worksize of c_matrix evaluation.
  pure elemental function c_matrix_worksize(this) result(res)
    type(c_matrix), intent(in) :: this
    !! this :: c_matrix
    integer(IK)                :: res
    res = this%nw
  end function c_matrix_worksize
!
!| Evaluation the C matrix; G matrix is also calculated at the same time.<br>
!  If nx>=ny C(cb,nx,ny), else C(cb,ny,nx)
  pure subroutine c_matrix_eval(this, b, ms, X, Y, C, W)
    type(c_matrix), intent(in)  :: this
    !! c_matrix
    type(mol_block), intent(in) :: b
    !! mol_block, b must match the one used for initialization.
    integer(IK), intent(in)     :: ms(*)
    !! mol symmetry array, associated with b.
    real(RK), intent(in)        :: X(*)
    !! reference coordinate
    real(RK), intent(in)        :: Y(*)
    !! target coordinate
    real(RK), intent(inout)     :: C(*)
    !! main memory
    real(RK), intent(inout)     :: W(*)
    !! work memory
    integer(IK)                 :: s, m, n, dm, px, gx, gy, wx, wy
!
    s = mol_block_nsym(b)
    m = mol_block_napm(b)
    n = mol_block_nmol(b)
!
    dm = mol_block_each_size(b)
    px = mol_block_pointer(b)
    gx = this%w
    gy = gx + n
    wx = gx
    wy = wx + dm
!
    call eval_g_matrix(dm, this%cb, n, X(px), Y(px), C(this%p), W(gx), W(gy))
!
    call eval_c_matrix(b, s, m, n, dm, this%cb, ms, X(px), Y(px), &
  &                    C(this%p), W(wx), W(wy))
!
  contains
!
!   trace of self correlation matrix
    pure subroutine eval_g_matrix(dm, cb, n, X, Y, C, GX, GY)
      integer(IK), intent(in)        :: dm, cb, n
      real(RK), intent(in)           :: X(dm, *), Y(dm, *)
      real(RK), intent(inout)        :: C(cb, *)
      real(RK), intent(inout)        :: GX(n), GY(n)
      integer(IK)                    :: i, j
!
      do concurrent(i=1:n)
        GX(i) = dot(dm, X(1, i), 1, X(1, i), 1)
      end do
      do concurrent(i=1:n)
        GY(i) = dot(dm, Y(1, i), 1, Y(1, i), 1)
      end do
!
      do concurrent(i=1:n, j=1:n)
        block
          integer(IK) :: ic
!         if nx>=ny, C(s,nx,ny). else C(s,ny,nx)
          ic = i + (j - 1) * n
          C(1, ic) = GX(i) + GY(j)
        end block
      end do
!
    end subroutine eval_g_matrix
!
!   get correlation matrix C = Y^t@X and optimal rotation R^t
    pure subroutine eval_c_matrix(b, s, m, n, dm, cb, ms, X, Y, C, WX, WY)
      type(mol_block), intent(in) :: b
      integer(IK), intent(in)     :: s, m, dm, n, cb, ms(*)
      real(RK), intent(in)        :: X(dm, *), Y(dm, *)
      real(RK), intent(inout)     :: C(cb, *)
      real(RK), intent(inout)     :: WX(dm), WY(dm, n)
      integer(IK)                 :: i, j
!
      call copy(dm * n, Y, 1, WY, 1)
!
      do i = 1, n
        call copy(dm, X(1, i), 1, WX, 1)
        do concurrent(j=1:n)
          block
            integer(IK) :: ic
            ic = i + (j - 1) * n
            call calc_cov(b, s, m, dm, ms, WX, WY(1, j), C(2, ic))
          end block
        end do
      end do
!
    end subroutine eval_c_matrix
!
    pure subroutine calc_cov(b, s, m, dm, ms, WX, WY, C)
      type(mol_block), intent(in) :: b
      integer(IK), intent(in)     :: s, m, dm, ms(*)
      real(RK), intent(in)        :: WX(D, *)
      real(RK), intent(inout)     :: WY(D, *)
      real(RK), intent(inout)     :: C(DD, *)
      integer(IK)                 :: i
!
        call gemm('N', 'T', D, D, m, ONE, WY, D, WX, D, ZERO, C(1, 1), D)
!
        do concurrent(i=2:s)
          call mol_block_swap(b, i - 1, ms, WY)
          call gemm('N', 'T', D, D, m, ONE, WY, D, WX, D, ZERO, C(1, i), D)
          call mol_block_inverse_swap(b, i - 1, ms, WY)
        end do
!
    end subroutine calc_cov
!
  end subroutine c_matrix_eval
!
!| Add CIJs to partial covariance matrix C.
  pure subroutine c_matrix_add(this, b, i, j, s, W, G, C)
    type(c_matrix), intent(in)  :: this
    !! this :: c_matrix
    type(mol_block), intent(in) :: b
    !! b    :: mol_block
    integer(IK), intent(in)     :: i
    !! i    :: row index
    integer(IK), intent(in)     :: j
    !! collumn index
    integer(IK), intent(in)     :: s
    !! symmetry index
    real(RK), intent(in)        :: W(*)
    !! main memory
    real(RK), intent(inout)     :: G
    !! partial auto variance matrix
    real(RK), intent(inout)     :: C(*)
    !! partial covariance matrix
    integer(IK)                 :: k
!
    k = this%p + this%cb * (i - 1) + this%cl * (j - 1)
    G = G + W(k)
    k = k + DD * (s - 1) + 1
    call axpy(DD, ONE, W(k), 1, C, 1)
!
  end subroutine c_matrix_add
!
end module mod_c_matrix

