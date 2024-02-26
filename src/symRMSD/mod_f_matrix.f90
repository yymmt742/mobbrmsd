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
  public :: memsize_f_matrix
  public :: worksize_f_matrix
  public :: f_matrix_eval
!
  !| f_matrix
  type f_matrix
    private
    sequence
    !| p  :: pointer to S.
    integer(IK), public :: p
    !| w  :: pointer to work.
    integer(IK), public :: w
    !| nn :: nx * ny
    integer(IK)         :: nn
    !| nw :: work memory size.
    integer(IK)         :: nw
  end type f_matrix
!
  interface f_matrix
    module procedure f_matrix_new
  end interface f_matrix
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
    res%nn = b%n1 * b%n2
    res%nw = worksize_sdmin()
!
  end function f_matrix_new
!
!| Inquire memsize of f_matrix.
  pure elemental function memsize_f_matrix(this) result(res)
    !| this :: f_matrix
    type(f_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%nn
  end function memsize_f_matrix
!
!| Inquire worksize of f_matrix.
  pure elemental function worksize_f_matrix(this) result(res)
    !| this :: f_matrix
    type(f_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%nw
  end function worksize_f_matrix
!
!| Evaluation the D matrix.<br>
!  If nx>=ny D(nx,ny), else D(ny,nx)
  pure subroutine f_matrix_eval(this, b, C, F, W)
    !| this :: f_matrix
    type(f_matrix), intent(in)  :: this
    !| b    :: mol_block
    type(mol_block), intent(in) :: b
    !| C    :: covariacne matrix C(cb, n1, n2)
    real(RK), intent(in)        :: C(*)
    !| F    :: main memory
    real(RK), intent(inout)     :: F(*)
    !| W    :: work array
    real(RK), intent(inout)     :: W(*)
    integer(IK)                 :: cb
!
    cb = 1 + DD * b%s
    call eval_f_matrix(cb, this%nn, C, F(this%p), W(this%w))
!
  contains
!
    pure subroutine eval_f_matrix(cb, nn, C, F, W)
      integer(IK), intent(in) :: cb, nn
      real(RK), intent(in)    :: C(cb, *)
      real(RK), intent(inout) :: F(*), W(*)
      integer(IK)             :: i, j
!
!!!   get squared displacement
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
  end subroutine f_matrix_eval
!
end module mod_f_matrix

