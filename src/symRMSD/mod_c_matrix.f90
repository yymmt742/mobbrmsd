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
    !| nn :: nx * ny
    integer(IK)         :: nn
    !| dm :: dm = D * m
    integer(IK)         :: dm
    !| cb :: number of elements in a sell. cb = DD * res%b%s + 1
    integer(IK)         :: cb
    !| cl :: number of elements in a line. cl = cb * MAX(nx, ny)
    integer(IK)         :: cl
    !| px :: pointer to x.
    integer(IK)         :: px
    !| py :: pointer to y.
    integer(IK)         :: py
    !| pc :: pointer to C.
    integer(IK)         :: pc
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
!| Constructer of c_matrix
  pure elemental function c_matrix_new(b, p) result(res)
    !| b :: mol_block
    type(mol_block), intent(in) :: b
    !| p :: pointer to C.
    integer(IK), intent(in)     :: p
    type(c_matrix)              :: res
!
    res%s = b%s
    res%m = b%m
    res%nx = b%x%n
    res%ny = b%x%n
    res%nn = res%nx * res%ny
    res%dm = D * res%m
    res%cb = 1 + DD * res%s
    res%cl = res%cb * MAX(res%nx, res%ny)
!
    res%px = b%x%p
    res%py = b%y%p
    res%pc = p
!
  end function c_matrix_new
!
!| Inquire memsize of c_matrix.
  pure elemental function c_matrix_memsize(this) result(res)
    !| this :: c_matrix
    type(c_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%cb * this%nn
  end function c_matrix_memsize
!
!| Evaluation C_matrix.
  pure subroutine c_matrix_eval(this, ms, W)
    !| this :: c_matrix
    type(c_matrix), intent(in)      :: this
    !| ms   :: mol_symmetry
    class(mol_symmetry), intent(in) :: ms
    !| W    :: work memory
    real(RK), intent(inout)         :: W(*)
!
    call eval(this%s, this%m, this%nx, this%ny, this%nn, this%dm, this%cb, ms, &
  &           W(this%px), W(this%py), W(this%pc))
!
  contains
!
    pure subroutine eval(s, m, nx, ny, nn, dm, cb, ms, X, Y, C)
      integer(IK), intent(in)        :: s, m, nx, ny, nn
      integer(IK), intent(in)        :: dm, cb
      type(mol_symmetry), intent(in) :: ms
      real(RK), intent(in)           :: X(dm, nx)
      real(RK), intent(in)           :: Y(dm, ny)
      real(RK), intent(inout)        :: C(cb, nn)
      real(RK)                       :: GX(nx), GY(ny)
      integer(IK)                    :: i, j, k
!
      do concurrent(j=1:nx)
        GX(j) = DDOT(dm, X(1, j), 1, X(1, j), 1)
      end do
!
      do concurrent(j=1:ny)
        GY(j) = DDOT(dm, Y(1, j), 1, Y(1, j), 1)
      end do
!
      do concurrent(i=0:s - 1, j=1:nx, k=1:ny)
        block
          integer(IK) :: is, ic
!         if nx>=ny, C(s,nx,ny). else C(s,ny,nx)
          is = 2 + i * DD
          ic = MERGE(2 + j + (k - 1) * nx, 2 + k + (j - 1) * ny, nx >= ny)
!         trace of self correlation matrix
          C(1, ic) = GX(j) + GY(k)
!         get correlation matrix C = Y^t@X and optimal rotation R^t
          call calc_cov(m, dm, i, ms, X(1, j), Y(1, k), C(is, ic))
        end block
      end do
!
    end subroutine eval
!
    pure subroutine calc_cov(m, dm, si, ms, X, Y, C)
      integer(IK), intent(in)        :: m, dm, si
      type(mol_symmetry), intent(in) :: ms
      real(RK), intent(in)           :: X(*), Y(*)
      real(RK), intent(inout)        :: C(*)
      real(RK)                       :: WX(dm), WY(dm)
!
        call DCOPY(dm, X, 1, WX, 1)
        call DCOPY(dm, Y, 1, WY, 1)
        call ms%swap(D, WY, si)
        call DGEMM('N', 'T', D, D, m, ONE, WY, D, WX, D, ZERO, C, D)
!
    end subroutine calc_cov
!
  end subroutine c_matrix_eval
!
end module mod_c_matrix

