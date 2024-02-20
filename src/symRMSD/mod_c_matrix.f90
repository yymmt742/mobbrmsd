!| Module for manage C matrix.<br>
!  {C_IJs} :: Covariance matrices.<br>
!    C_IJs(d,d) = Y_J @ Q_s @ X_I^T<br>
!    - X_I :: I-th molecule in X.<br>
!    - Y_J :: J-th molecule in Y.<br>
!    - Q_s :: Permutation matrix on m.
module mod_c_matrix
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_mol_symmetry
  use mod_mol_block
  implicit none
  private
  public :: c_matrix
  public :: c_matrix_memsize
  public :: c_matrix_worksize
  public :: c_matrix_init
  public :: c_matrix_eval
!
  !| c_matrix<br>
  !  C_IJs is the d x d matrix and {C_IJs} is the third-order tensor of nx x ny x s.<br>
  !  To quickly find the rotation matrix, C is stored with G_IJ given by<br>
  !    G_IJ = Tr[X_I @ X_I^T] + Tr[Y_J @ Y_J^T].<br>
  !  G does not change with respect to s.<br>
  !  If nx >= ny, {C_IJs} is stored as C(cb,nx,ny), otherwise, C(cb,ny,nx), <br>
  !  where cb = 1 + s * d * d.<br>
  !  C(:,I,J) = [G_IJ, C_IJ1, C_IJ2, ..., C_IJS] with C_IJs(D,D) := Y_J @ Q_s @ X_I^T.
  type c_matrix
    private
    sequence
    !| s  :: number of molecular symmetry.
    integer(IK)         :: s
    !| m  :: number of atoms in a molecule.
    integer(IK)         :: m
    !| nx :: number of molecule in X.
    integer(IK)         :: nx
    !| ny :: number of molecule in Y.
    integer(IK)         :: ny
    !| cb :: number of elements in a sell. cb = DD * res%b%s + 1
    integer(IK), public :: cb
    !| cl :: number of elements in a line. cl = cb * MAX(nx, ny)
    integer(IK), public :: cl
    !| nl :: number of row. nl = MIN(nx, ny)
    integer(IK), public :: nl
    !| px :: pointer to x.
    integer(IK)         :: px
    !| py :: pointer to y.
    integer(IK)         :: py
    !| gx :: pointer to GX.
    integer(IK)         :: gx
    !| gy :: pointer to GY.
    integer(IK)         :: gy
    !| wx :: pointer to WX.
    integer(IK)         :: wx
    !| wy :: pointer to WY.
    integer(IK)         :: wy
    !| pc :: pointer to C.
    integer(IK), public :: pc
  end type c_matrix
!
  interface c_matrix
    module procedure c_matrix_new
  end interface c_matrix
!
  interface
    include 'dgemm.h'
    include 'ddot.h'
    include 'dcopy.h'
  end interface
!
contains
!
!| Constructer
  pure elemental function c_matrix_new(b) result(res)
    !| b :: mol_block, must be initialized.
    type(mol_block), intent(in) :: b
    type(c_matrix)              :: res
!
    res%s = b%s
    res%m = b%m
    res%nx = b%x%n
    res%ny = b%x%n
    res%cb = 1 + DD * res%s
    res%cl = res%cb * MAX(res%nx, res%ny)
    res%nl = MIN(res%nx, res%ny)
!
    res%px = b%x%p
    res%py = b%y%p
!
    call c_matrix_init(res)
!
  end function c_matrix_new
!
!| Inquire memsize of c_matrix.
  pure elemental function c_matrix_memsize(this) result(res)
    !| this :: c_matrix
    type(c_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%cb * this%nx * this%ny
  end function c_matrix_memsize
!
!| Inquire worksize of c_matrix.
  pure elemental function c_matrix_worksize(this) result(res)
    !| this :: c_matrix
    type(c_matrix), intent(in) :: this
    integer(IK)                :: res
    res = MAX(D * this%m * (1 + this%ny), this%nx + this%ny)
  end function c_matrix_worksize
!
!| Initializer of c_matrix array.<br>
!  By default, memory is allocated as follows.<br>
!  |-C1-|----W1----|<br>
!  |----|--C2--|--------W2-------|<br>
!  |-----------|--C3--|-W3-|<br>
!  Therefore, the maximum memory allocation size is MAX( SUM_i^I |Ci| + |W_I| ).
  pure subroutine c_matrix_init(this, p)
    type(c_matrix), intent(inout)     :: this
    integer(IK), intent(in), optional :: p
    integer(IK)                       :: q
!
    if (PRESENT(p)) then; q = p
    else; q = 1
    end if
!
    this%pc = q
    q = q + c_matrix_memsize(this)
    this%gx = q
    this%gy = this%gx + this%nx
    this%wx = q
    this%wy = this%wx + D * this%m
!
  end subroutine c_matrix_init
!
!| Evaluation the C matrix; G matrix is also calculated at the same time.<br>
!  If nx>=ny C(cb,nx,ny), else C(cb,ny,nx)
  pure subroutine c_matrix_eval(this, ms, X, Y, W)
    !| this :: c_matrix
    type(c_matrix), intent(in)      :: this
    !| ms   :: mol_symmetry
    class(mol_symmetry), intent(in) :: ms
    !| X    :: reference coordinate
    real(RK), intent(inout)         :: X(*)
    !| Y    :: target coordinate
    real(RK), intent(inout)         :: Y(*)
    !| W    :: work memory
    real(RK), intent(inout)         :: W(*)
    integer(IK)                     :: dm
!
    dm = D * this%m
!
    call eval_g_matrix(dm, this%cb, this%nx, this%ny, &
   &                   X(this%px), Y(this%py), W(this%pc), W(this%gx), W(this%gy))
!
    call eval_c_matrix(this%s, this%m, dm, this%nx, this%ny, this%cb, ms, &
  &                    X(this%px), Y(this%py), W(this%pc), W(this%wx), W(this%wy))
!
  contains
!
!   trace of self correlation matrix
    pure subroutine eval_g_matrix(dm, cb, nx, ny, X, Y, C, GX, GY)
      integer(IK), intent(in)        :: dm, cb, nx, ny
      real(RK), intent(in)           :: X(dm, *), Y(dm, *)
      real(RK), intent(inout)        :: C(cb, *)
      real(RK), intent(inout)        :: GX(nx), GY(ny)
      integer(IK)                    :: i, j
!
      do concurrent(i=1:nx)
        GX(i) = DDOT(dm, X(1, i), 1, X(1, i), 1)
      end do
      do concurrent(i=1:ny)
        GY(i) = DDOT(dm, Y(1, i), 1, Y(1, i), 1)
      end do
!
      do concurrent(i=1:nx, j=1:ny)
        block
          integer(IK) :: ic
!         if nx>=ny, C(s,nx,ny). else C(s,ny,nx)
          ic = MERGE(i + (j - 1) * nx, j + (i - 1) * ny, nx >= ny)
          C(1, ic) = GX(i) + GY(j)
        end block
      end do
!
    end subroutine eval_g_matrix
!
!   get correlation matrix C = Y^t@X and optimal rotation R^t
    pure subroutine eval_c_matrix(s, m, dm, nx, ny, cb, ms, X, Y, C, WX, WY)
      integer(IK), intent(in)        :: s, m, dm, nx, ny, cb
      type(mol_symmetry), intent(in) :: ms
      real(RK), intent(in)           :: X(dm, nx), Y(dm, ny)
      real(RK), intent(inout)        :: C(cb, *)
      real(RK), intent(inout)        :: WX(dm), WY(dm, ny)
      integer(IK)                    :: j, k
!
      call DCOPY(dm * ny, Y, 1, WY, 1)
!
      do j = 1, nx
        call DCOPY(dm, X(1, j), 1, WX, 1)
        do concurrent(k=1:ny)
          block
            integer(IK) :: ic
            ic = MERGE(j + (k - 1) * nx, k + (j - 1) * ny, nx >= ny)
            call calc_cov(s, m, dm, ms, WX, WY(1, k), C(2, ic))
          end block
        end do
      end do
!
    end subroutine eval_c_matrix
!
    pure subroutine calc_cov(s, m, dm, ms, WX, WY, C)
      integer(IK), intent(in)        :: s, m, dm
      type(mol_symmetry), intent(in) :: ms
      real(RK), intent(in)           :: WX(*)
      real(RK), intent(inout)        :: WY(*)
      real(RK), intent(inout)        :: C(DD, *)
      integer(IK)                    :: i
!
        call DGEMM('N', 'T', D, D, m, ONE, WY, D, WX, D, ZERO, C(1, 1), D)
!
        do concurrent(i=2:s)
          call ms%swap(D, WY, i - 1)
          call DGEMM('N', 'T', D, D, m, ONE, WY, D, WX, D, ZERO, C(1, i), D)
          call ms%reverse(D, WY, i - 1)
        end do
!
    end subroutine calc_cov
!
  end subroutine c_matrix_eval
!
end module mod_c_matrix

