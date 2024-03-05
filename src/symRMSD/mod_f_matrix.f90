!| Module for manage S matrix.<br>
!  S(nx, ny) :: Surplus matrices.<br>
!    - S_IJ = min_{R,s} Tr[C_IJs @ R]<br>
module mod_f_matrix
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_c_matrix
  use mod_mol_block
  use mod_rotation_matrix
  use mod_Hungarian
  implicit none
  private
  public :: f_matrix
  public :: f_matrix_tuple
  public :: f_matrix_memsize
  public :: f_matrix_worksize
  public :: f_matrix_eval
!
!| f_matrix <br>
!  - F is the n x n matrix. <br>
!  - F(I,J) = min_{R, s} tr[R@C_IJS] <br>
  type f_matrix
    private
    sequence
    !| p  :: pointer to S.
    integer(IK), public :: p
    !| w  :: pointer to work.
    integer(IK), public :: w
    !| nn :: n * n
    integer(IK)         :: nn
    !| nw :: work memory size.
    integer(IK)         :: nw
  end type f_matrix
!
  interface f_matrix
    module procedure f_matrix_new
  end interface f_matrix
!
!| A set of f_matrix and work arrays. <br>
!  This is mainly used for passing during initialization.
  type f_matrix_tuple
    type(f_matrix)        :: f
    !! header
    real(RK), allocatable :: x(:)
    !! main memory.
    real(RK), allocatable :: w(:)
    !! work memory.
  contains
    final :: f_matrix_tuple_destroy
  end type f_matrix_tuple
!
  interface f_matrix_tuple
    module procedure f_matrix_tuple_new
  end interface f_matrix_tuple
!
contains
!
!| Constructer
  pure elemental function f_matrix_new(b) result(res)
    !| b :: mol_block, must be initialized.
    type(mol_block), intent(in) :: b
    type(f_matrix)              :: res
!
    res%p  = 1
    res%w  = 1
    res%nn = mol_block_nmol(b)**2
    res%nw = sdmin_worksize()
!
  end function f_matrix_new
!
!| Constructer
  pure elemental function f_matrix_tuple_new(b) result(res)
    type(mol_block), intent(in) :: b
    !! mol_block, must be initialized.
    type(f_matrix_tuple)        :: res
!
    res%f = f_matrix(b)
    allocate (res%x(f_matrix_memsize(res%f)))
    allocate (res%w(f_matrix_worksize(res%f)))
!
  end function f_matrix_tuple_new
!
!| Inquire memsize of f_matrix.
  pure elemental function f_matrix_memsize(this) result(res)
    !| this :: f_matrix
    type(f_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%nn
  end function f_matrix_memsize
!
!| Inquire worksize of f_matrix.
  pure elemental function f_matrix_worksize(this) result(res)
    !| this :: f_matrix
    type(f_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%nw
  end function f_matrix_worksize
!
!| Evaluation the D matrix.<br>
!  If nx>=ny D(nx,ny), else D(ny,nx)
  pure subroutine f_matrix_eval(this, b, c, X, F, W)
    type(f_matrix), intent(in)  :: this
    !! f_matrix
    type(mol_block), intent(in) :: b
    !! mol_block, b must match the one used for initialization.
    type(c_matrix), intent(in)  :: c
    !! covariacne matrix c.
    real(RK), intent(inout)     :: X(*)
    !! main memory of C.
    real(RK), intent(inout)     :: F(*)
    !! main memory of F.
    real(RK), intent(inout)     :: W(*)
    !! work array.
    integer(IK)                 :: cb
!
    cb = c_matrix_blocksize(c)
    call eval_f_matrix(cb, this%nn, X(c%p), F(this%p), W(this%w))
!
  end subroutine f_matrix_eval
!
  pure subroutine eval_f_matrix(cb, nn, C, F, W)
    integer(IK), intent(in) :: cb, nn
    real(RK), intent(in)    :: C(cb, *)
    real(RK), intent(inout) :: F(*), W(*)
    integer(IK)             :: i, j
!
!   get squared displacement
    do j = 1, nn
      F(j) = RHUGE
      do i = 2, cb, DD
        call estimate_sdmin(C(1, j), C(i, j), W)
        F(j) = MIN(F(j), W(1))
      end do
    end do
!
  end subroutine eval_f_matrix
!
  pure elemental subroutine f_matrix_tuple_destroy(this)
    type(f_matrix_tuple), intent(inout) :: this
    if (ALLOCATED(this%x)) deallocate (this%x)
    if (ALLOCATED(this%w)) deallocate (this%w)
  end subroutine f_matrix_tuple_destroy
!
end module mod_f_matrix

