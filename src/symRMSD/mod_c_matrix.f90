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
  public :: c_matrix_memsize
  public :: c_matrix_worksize
  public :: c_matrix_blocksize
  public :: c_matrix_eval
  public :: c_matrix_add
!
  integer(IK), parameter :: header_size = 4
!
  integer(IK), parameter :: nl = 1
  !! number of row. nl = n.
  integer(IK), parameter :: cb = 2
  !! number of elements in a sell. cb = DD * res%b%s + 1.
  integer(IK), parameter :: cl = 3
  !! number of elements in a line. cl = cb * n.
  integer(IK), parameter :: nw = 4
  !! number of work array. nw = MAX(dmn + dm, n + n)
!
!| c_matrix <br>
!  - C_IJs is the d x d matrix and {C_IJs} is the third-order tensor of nx x ny x s. <br>
!  - To quickly find the rotation matrix, C is stored with G_IJ given by <br>
!  --- G_IJ = Tr[X_I @ X_I^T] + Tr[Y_J @ Y_J^T]. <br>
!  - G does not change with respect to s. <br>
!  - C(:,I,J) = [G_IJ, C_IJ1, C_IJ2, ..., C_IJS] with C_IJs(D,D) := Y_J @ Q_s @ X_I^T. <br>
!  This is mainly used for passing during initialization.
  type c_matrix
    integer(IK)              :: q(header_size)
    !! header
  end type c_matrix
!
  interface c_matrix
    module procedure c_matrix_new
  end interface c_matrix
!
contains
!
!| Constructer
  pure function c_matrix_new(b) result(res)
    integer(IK), intent(in) :: b(*)
    !! mol_block, must be initialized.
    type(c_matrix)    :: res
!
    res%q(cb) = 1 + DD * mol_block_nsym(b)
    res%q(nl) = mol_block_nmol(b)
    res%q(cl) = res%q(cb) * res%q(nl)
    res%q(nw) = MAX((res%q(nl) + 1) * mol_block_each_size(b), res%q(nl) * 2)
!
  end function c_matrix_new
!
!| Inquire blocksize of c_matrix.
  pure function c_matrix_blocksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! c_matrix
    integer(IK)             :: res
    res = q(cb)
  end function c_matrix_blocksize
!
!| Inquire memsize of c_matrix.
  pure function c_matrix_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! c_matrix
    integer(IK)             :: res
    res = q(cl) * q(nl)
  end function c_matrix_memsize
!
!| Inquire worksize of c_matrix evaluation.
  pure function c_matrix_worksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! c_matrix
    integer(IK)             :: res
    res = q(nw)
  end function c_matrix_worksize
!
!| Evaluation the C matrix; G matrix is also calculated at the same time.<br>
!  If nx>=ny C(cb,nx,ny), else C(cb,ny,nx)
  pure subroutine c_matrix_eval(q, b, X, Y, C, W)
    integer(IK), intent(in) :: q(*)
    !! c_matrix
    integer(IK), intent(in) :: b(*)
    !! mol symmetry array, associated with b.
    real(RK), intent(in)    :: X(*)
    !! reference coordinate
    real(RK), intent(in)    :: Y(*)
    !! target coordinate
    real(RK), intent(inout) :: C(*)
    !! main memory
    real(RK), intent(inout) :: W(*)
    !! work memory
    integer(IK), parameter  :: gx = 1
    integer(IK), parameter  :: wx = 1
    integer(IK)             :: s, m, n, dm, gy, wy
!
    s = mol_block_nsym(b)
    m = mol_block_napm(b)
    n = mol_block_nmol(b)
!
    dm = mol_block_each_size(b)
    gy = gx + n
    wy = wx + dm
!
    call eval_g_matrix(dm, q(cb), n, X, Y, C, W(gx), W(gy))
    call eval_c_matrix(b, s, m, n, dm, q(cb), X, Y, C, W(wx), W(wy))
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
    pure subroutine eval_c_matrix(b, s, m, n, dm, cb, X, Y, C, WX, WY)
      integer(IK), intent(in)     :: b(*), s, m, dm, n, cb
      real(RK), intent(in)        :: X(dm, *), Y(dm, *)
      real(RK), intent(inout)     :: C(cb, n, *)
      real(RK), intent(inout)     :: WX(dm), WY(dm, n)
      integer(IK)                 :: i, j
!
      call copy(dm * n, Y, 1, WY, 1)
!
      do j = 1, n
        call copy(dm, X(1, j), 1, WX, 1)
        do concurrent(i=1:n)
          call calc_cov(b, s, m, dm, WX, WY(1, i), C(2, i, j))
        end do
      end do
!
    end subroutine eval_c_matrix
!
    pure subroutine calc_cov(b, s, m, dm, WX, WY, C)
      integer(IK), intent(in)     :: b(*), s, m, dm
      real(RK), intent(in)        :: WX(D, m)
      real(RK), intent(inout)     :: WY(D, m)
      real(RK), intent(inout)     :: C(DD, *)
      integer(IK)                 :: i
!
        call gemm('N', 'T', D, D, m, ONE, WY, D, WX, D, ZERO, C(1, 1), D)
!
        do concurrent(i=2:s)
          call mol_block_swap(b, i - 1, WY)
          call gemm('N', 'T', D, D, m, ONE, WY, D, WX, D, ZERO, C(1, i), D)
          call mol_block_inverse_swap(b, i - 1, WY)
        end do
!
    end subroutine calc_cov
!
  end subroutine c_matrix_eval
!
!| Add CIJs to partial covariance matrix C.
  pure subroutine c_matrix_add(q, i, j, s, C, G, Cp)
    integer(IK), intent(in) :: q(*)
    !! this :: c_matrix
    integer(IK), intent(in) :: i
    !! i    :: row index
    integer(IK), intent(in) :: j
    !! collumn index
    integer(IK), intent(in) :: s
    !! symmetry index
    real(RK), intent(in)    :: C(*)
    !! main memory
    real(RK), intent(inout) :: G
    !! partial auto variance matrix
    real(RK), intent(inout) :: Cp(*)
    !! partial covariance matrix
    integer(IK)             :: k
!
    k = q(cl) * (i - 1) + q(cb) * (j - 1) + 1
    G = G + C(k)
    k = k + 1 + DD * (s - 1)
    call axpy(DD, ONE, C(k), 1, Cp, 1)
!
  end subroutine c_matrix_add
!
end module mod_c_matrix

