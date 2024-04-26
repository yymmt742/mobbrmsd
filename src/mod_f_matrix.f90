!| Module for managing f-matrix, free rotation cost matrix, F(M, M). <br>
!  \[ \mathbf{F}_{IJ} = \min_{\mathbf{R},s}\text{Tr}\left[\mathbf{C}_{IJs} \mathbf{R} \right] \] <br>
!  \( \mathbf{C}_{IJs} \) :: Covariance matrix of \( \mathbf{X}_I \) and \( \mathbf{Y}_I \) with \( s \)-th molecular permutation.<br>
!  \( \mathbf{R} \) :: Rotation matrix on \( \mathbb{R}^{d\times d} \).<br>
module mod_f_matrix
  use blas_lapack_interface, only: D, DD
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_c_matrix
  use mod_mol_block
  use mod_rotation
  use mod_Hungarian
  implicit none
  private
  public :: f_matrix
  public :: f_matrix_memsize
  public :: f_matrix_worksize
  public :: f_matrix_eval
!
  integer(IK), parameter :: header_size = 1
  integer(IK), parameter :: nn = 1
!
!| f_matrix <br>
!  This is mainly used for passing during initialization.
  type f_matrix
    integer(IK)           :: q(header_size)
    !! header
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
    !! mol_block, must be initialized
    type(f_matrix)          :: res
!
    res%q(nn) = mol_block_nmol(b)**2
!
  end function f_matrix_new
!
!| Inquire memsize of f_matrix
  pure function f_matrix_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! f_matrix
    integer(IK)             :: res
    res = q(nn)
  end function f_matrix_memsize
!
!| Inquire worksize of f_matrix
  pure function f_matrix_worksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! f_matrix
    integer(IK)             :: res
    res = q(nn) * (sdmin_worksize() + 1)
  end function f_matrix_worksize
!
!| Evaluation the D matrix
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
    integer(IK)             :: i, cb, nw
!
    cb = c_matrix_blocksize(qc)
    nw = sdmin_worksize() + 1
!
    do concurrent(i=1:q(nn))
      block
        integer(IK) :: ic, iw
        ic = cb * (i - 1) + 1
        iw = nw * (i - 1) + 1
        call eval_f_matrix(cb, C(ic), W(iw))
        F(i) = W(iw)
      end block
    end do
!
  end subroutine f_matrix_eval
!
  pure subroutine eval_f_matrix(cb, C, W)
    integer(IK), intent(in) :: cb
    real(RK), intent(in)    :: C(cb)
    real(RK), intent(inout) :: W(*)
    integer(IK)             :: i
!   get - max tr[CR]
    W(1) = RHUGE
    do i = 2, cb, DD
      call estimate_rcmax(C(1), C(i), W(2))
      W(1) = MIN(W(1), -W(2))
    end do
  end subroutine eval_f_matrix
!
end module mod_f_matrix

