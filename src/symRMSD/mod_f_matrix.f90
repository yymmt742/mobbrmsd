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
  public :: f_matrix_memsize
  public :: f_matrix_worksize
  public :: f_matrix_eval
!
  integer(IK), parameter :: header_size = 2
  integer(IK), parameter :: nn = 1
  !! n * n
  integer(IK), parameter :: nw = 2
  !! work memory size.
!
!| f_matrix <br>
!  - F is the n x n matrix. <br>
!  - F(I,J) = min_{R, s} tr[R@C_IJS] <br>
!  This is mainly used for passing during initialization.
  type f_matrix
    integer(IK)           :: q(header_size)
    !! header
    real(RK), allocatable :: x(:)
    !! main memory.
    real(RK), allocatable :: w(:)
    !! work memory.
  contains
    final :: f_matrix_destroy
  end type f_matrix
!
  interface f_matrix
    module procedure f_matrix_new
  end interface f_matrix
!
contains
!
!| Constructer
  pure function f_matrix_new(b) result(res)
    integer(IK), intent(in) :: b(*)
    !! mol_block, must be initialized.
    type(f_matrix)    :: res
!
    res%q(nn) = mol_block_nmol(b)**2
    res%q(nw) = sdmin_worksize()

    allocate (res%x(f_matrix_memsize(res%q)))
    allocate (res%w(f_matrix_worksize(res%q)))
!
  end function f_matrix_new
!
!| Inquire memsize of f_matrix.
  pure function f_matrix_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! f_matrix
    integer(IK)             :: res
    res = q(nn)
  end function f_matrix_memsize
!
!| Inquire worksize of f_matrix.
  pure function f_matrix_worksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! f_matrix
    integer(IK)             :: res
    res = q(nw)
  end function f_matrix_worksize
!
!| Evaluation the D matrix.<br>
!  If nx>=ny D(nx,ny), else D(ny,nx)
  pure subroutine f_matrix_eval(q, qc, C, F, W)
    integer(IK), intent(in) :: q(*)
    !! header of f_matrix
    integer(IK), intent(in) :: qc(*)
    !! header of covariacne matrix C.
    real(RK), intent(inout) :: C(*)
    !! covariacne matrix c.
    real(RK), intent(inout) :: F(*)
    !! main memory of F.
    real(RK), intent(inout) :: W(*)
    !! work array.
    integer(IK)             :: cb
!
    cb = c_matrix_blocksize(qc)
    call eval_f_matrix(cb, q(nn), C, F, W)
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
  pure elemental subroutine f_matrix_destroy(this)
    type(f_matrix), intent(inout) :: this
    if (ALLOCATED(this%x)) deallocate (this%x)
    if (ALLOCATED(this%w)) deallocate (this%w)
  end subroutine f_matrix_destroy
!
end module mod_f_matrix

