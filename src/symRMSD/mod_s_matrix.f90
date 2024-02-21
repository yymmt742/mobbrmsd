!| Module for manage S matrix.<br>
!  S(nx, ny) :: Surplus matrices.<br>
!    - S_IJ = min_{R,s} Tr[C_IJs @ R]<br>
module mod_s_matrix
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_c_matrix
  use mod_mol_block
  use mod_rotation_matrix
  implicit none
  private
  public :: s_matrix
  public :: memsize_s_matrix
  public :: worksize_s_matrix
  public :: s_matrix_eval
!
  !| s_matrix
  type s_matrix
    private
    sequence
    !| p  :: pointer to S.
    integer(IK), public :: p
    !| w  :: pointer to work.
    integer(IK), public :: w
    !| nn :: memory size.
    integer(IK)         :: nn
    !| nw :: work memory size.
    integer(IK)         :: nw
  end type s_matrix
!
  interface s_matrix
    module procedure s_matrix_new
  end interface s_matrix
!
contains
!
!| Constructer
  pure elemental function s_matrix_new(b) result(res)
    !| b :: mol_block, must be initialized.
    type(mol_block), intent(in) :: b
    type(s_matrix)              :: res
!
    res%p  = 1
    res%w  = 1
    res%nn = b%x%n * b%y%n
    res%nw = worksize_sdmin()
!
  end function s_matrix_new
!
!| Inquire memsize of s_matrix.
  pure elemental function memsize_s_matrix(this) result(res)
    !| this :: s_matrix
    type(s_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%nn
  end function memsize_s_matrix
!
!| Inquire worksize of s_matrix.
  pure elemental function worksize_s_matrix(this) result(res)
    !| this :: s_matrix
    type(s_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%nw
  end function worksize_s_matrix
!
!| Evaluation the D matrix.<br>
!  If nx>=ny D(nx,ny), else D(ny,nx)
  pure subroutine s_matrix_eval(this, b, C, S, W)
    !| this :: s_matrix
    type(s_matrix), intent(in)  :: this
    !| b    :: mol_block
    type(mol_block), intent(in) :: b
    !| C    :: covariacne matrix C(cb, n1, n2)
    real(RK), intent(in)        :: C(*)
    !| S    :: main memory
    real(RK), intent(inout)     :: S(*)
    !| W    :: work array
    real(RK), intent(inout)     :: W(*)
    integer(IK)                 :: cb
!
    cb = 1 + DD * b%s
    call eval_s_matrix(cb, this%nn, C, S(this%p), W(this%w))
!
  contains
!
    pure subroutine eval_s_matrix(cb, nn, C, S, W)
      integer(IK), intent(in) :: cb, nn
      real(RK), intent(in)    :: C(cb, *)
      real(RK), intent(inout) :: S(*), W(*)
      integer(IK)             :: i, j
!
!!!   get squared displacement
      do j = 1, nn
        S(j) = RHUGE
        do i = 2, cb, DD
          call estimate_sdmin(C(1, j), C(i, j), W)
          S(j) = MIN(S(j), W(1))
        end do
      end do
!
    end subroutine eval_s_matrix
!
  end subroutine s_matrix_eval
!
end module mod_s_matrix

