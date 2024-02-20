!| Module for manage D matrix.<br>
!  D(nx, ny) :: Residue matrices.<br>
!    - D_IJ = min_{R,s} Tr[C_IJs @ R]<br>
module mod_d_matrix
  use mod_params, only: D, DD, IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_c_matrix
  use mod_mol_block
  use mod_rotation_matrix
  implicit none
  private
  public :: d_matrix
  public :: memsize_d_matrix
  public :: worksize_d_matrix
  public :: d_matrix_eval
!
  !| d_matrix
  type d_matrix
    private
    sequence
    !| pd :: pointer to D.
    integer(IK), public :: p
    !| nm :: memory size
    integer(IK)         :: nn
    !| pd :: work memory size
    integer(IK)         :: nw
  end type d_matrix
!
  interface d_matrix
    module procedure d_matrix_new
  end interface d_matrix
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
  pure elemental function d_matrix_new(b) result(res)
    !| b :: mol_block, must be initialized.
    type(mol_block), intent(in) :: b
    type(d_matrix)              :: res
!
    res%p  = 1
    res%nn = b%x%n * b%y%n
    res%nw = worksize_sdmin() + 1
!
  end function d_matrix_new
!
!| Inquire memsize of d_matrix.
  pure elemental function memsize_d_matrix(this) result(res)
    !| this :: d_matrix
    type(d_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%nn
  end function memsize_d_matrix
!
!| Inquire worksize of d_matrix.
  pure elemental function worksize_d_matrix(this) result(res)
    !| this :: d_matrix
    type(d_matrix), intent(in) :: this
    integer(IK)                :: res
    res = this%nw
  end function worksize_d_matrix
!
!| Evaluation the D matrix.<br>
!  If nx>=ny D(nx,ny), else D(ny,nx)
  pure subroutine d_matrix_eval(this, b, cm, W)
    !| this :: d_matrix
    type(d_matrix), intent(in)      :: this
    !| b    :: mol_block
    type(mol_block), intent(in)     :: b
    !| cm   :: c_matrix
    type(c_matrix), intent(in)      :: cm
    real(RK), intent(inout)         :: W(*)
!
    call eval_d_matrix(b%s, cm%cb, this%nn, W(cm%p), W(this%p))
!
  contains
!
    pure subroutine eval_d_matrix(s, cb, nn, C, D)
      integer(IK), intent(in) :: s, cb, nn
      real(RK), intent(in)    :: C(cb, *)
      real(RK), intent(inout) :: D(*)
      integer(IK)             :: i, j
!
!!!   get squared displacement
      do i = 1, nn
        D(i) = RHUGE
        do j = 2, cb, DD
          call estimate_sdmin(C(1, i), C(j, i), D(i + 1))
          D(i) = MIN(D(i), D(i + 1))
        end do
      end do
!
    end subroutine eval_d_matrix
!
  end subroutine d_matrix_eval
!
end module mod_d_matrix

