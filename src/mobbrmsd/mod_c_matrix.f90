!| A module for managing c-matrices, tensor of covariance matrices. <br>
!  \( \{\mathbf{C}_{IJs}\}_{IJs} \) is third-order ( \(M\times M\times s \) ) tensor of matrix \( \mathbf{C}_{IJs}\in\mathbb{R}^{d\times d} \),
!  defined by<br>
!  \[ \mathbf{C}_{IJs} = \mathbf{Y}_J \mathbf{Q}_s \mathbf{X}_I^\top \] <br>
!  \( \mathbf{X}_I \) :: \( I \)-th molecule in reference coordinate, \( \mathbf{X} \in\mathbb{R}^{d\times n}\).<br>
!  \( \mathbf{Y}_J \) :: \( J \)-th molecule in target coordinate,    \( \mathbf{Y} \in\mathbb{R}^{d\times n} \).<br>
!  \( \mathbf{Q}_s \) :: Molecular permutation matrix on \( n \). <br>
!  To quickly find the rotation matrix, \( \mathbf{C}_{IJs} \) is stored with autocorrelation \( G_{IJ} \), defined by <br>
!  \[ G_{IJ} = \text{Tr}\left[\mathbf{X}_I\mathbf{X}_I^\top\right] + \text{Tr}\left[\mathbf{Y}_J\mathbf{Y}_J^\top\right] \] <br>
!  @note
!    \( G_{IJ} \) does not change with respect to molecular symmetry permutation index, \( s \). <br>
!    Therefore, data blocks are stored in three-dimensional array C(bs,M,M)
!    with the leading dimension with \( \text{bs}=1 + Sd^2 \). <br>
!    Data blocks are defined by \( \left[ G_{IJ}, \mathbf{C}_{IJ1}, \mathbf{C}_{IJ2}, \dots, \mathbf{C}_{IJS} \right] \) <br>
!  @endnote
module mod_c_matrix
  use blas_lapack_interface, only: D, DD
  use mod_cov, only: covdot, covcopy
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO, RHUGE
  use mod_mol_block
  implicit none
  private
  public :: c_matrix
  public :: c_matrix_memsize
  public :: c_matrix_worksize
  public :: c_matrix_blocksize
  public :: c_matrix_autocorr
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
!| C matrix manager.<br>
!   @note
!   This type is mainly used for passing during initialization.
!   @endnote
  type c_matrix
    integer(IK) :: q(header_size)
    !! header
  end type c_matrix
!
!| Constructer
  interface c_matrix
    module procedure c_matrix_new
  end interface c_matrix
!
  interface
    include 'dgemm.h'
    include 'sgemm.h'
  end interface
!
contains
!
!| Constructer
  pure function c_matrix_new(b) result(res)
    integer(IK), intent(in) :: b(*)
    !! mol_block header, must be initialized.
    type(c_matrix)          :: res
    res%q(cb) = 1 + DD * mol_block_nsym(b)
    res%q(nl) = mol_block_nmol(b)
    res%q(cl) = res%q(cb) * res%q(nl)
    res%q(nw) = MERGE(MAX((res%q(nl) + 1) * mol_block_each_size(b), res%q(nl) * 2), 0, res%q(nl) > 0)
  end function c_matrix_new
!
!| Inquire blocksize of c_matrix.
  pure function c_matrix_blocksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! c_matrix header.
    integer(IK)             :: res
    res = q(cb)
  end function c_matrix_blocksize
!
!| Inquire memsize of c_matrix.
  pure function c_matrix_memsize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! c_matrix header.
    integer(IK)             :: res
    res = q(cl) * q(nl)
  end function c_matrix_memsize
!
!| Inquire worksize of c_matrix evaluation.
  pure function c_matrix_worksize(q) result(res)
    integer(IK), intent(in) :: q(*)
    !! c_matrix header.
    integer(IK)             :: res
    res = q(nw)
  end function c_matrix_worksize
!
!| Evaluation the c-matrix.<br>
!  g-matrix is also calculated at the same time.<br>
!  At the end of the calculation, save the c-matrix C(cb,M,M) to C(\*). <br>
!  workarray W(\*) must be larger than c_matrix_worksize(q). <br>
  pure subroutine c_matrix_eval(q, b, X, Y, CX, XY, C, W)
    integer(IK), intent(in) :: q(*)
    !! c_matrix header.
    integer(IK), intent(in) :: b(*)
    !! mol block header.
    real(RK), intent(in)    :: X(*)
    !! reference coordinate, \(\mathbf{X}\). Arrays must be stored as X(D, n_apm, n_mol).
    real(RK), intent(in)    :: Y(*)
    !! target coordinate, \(\mathbf{Y}\). Arrays must be stored as Y(D, n_apm, n_mol).
    real(RK), intent(in)    :: CX(*)
    !! centroid of \(\mathbf{X}\).
    real(RK), intent(in)    :: CY(*)
    !! centroid of \(\mathbf{Y}\).
    real(RK), intent(inout) :: C(*)
    !! main memory
    real(RK), intent(inout) :: W(*)
    !! work memory
    integer(IK), parameter  :: gx = 1
    integer(IK), parameter  :: wy = 1
    integer(IK)             :: s, m, n, dm, gy, wx
!
    s = mol_block_nsym(b)
    m = mol_block_napm(b)
    n = mol_block_nmol(b)
!
    if (m < 1) then
      call zfill(q(cl) * q(nl), C, 1)
      return
    end if
!
    dm = mol_block_each_size(b)
    gy = gx + n
    wx = wy + dm
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
        GX(i) = dot(dm, X(1, i), X(1, i))
      end do
      do concurrent(i=1:n)
        GY(i) = dot(dm, Y(1, i), Y(1, i))
      end do
!
      do concurrent(i=1:n, j=1:n)
        block
          integer(IK) :: ic
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
      real(RK), intent(inout)     :: C(cb, n, *), WX(dm, n), WY(dm)
      integer(IK)                 :: i
!
      call copy(dm * n, X, WX)
!
      do i = 1, n
        call copy(dm, Y(1, i), WY)
        call calc_covariance(b, s, m, n, cb, WX, WY, C(1, 1, i))
      end do
!
    end subroutine eval_c_matrix
!
    pure subroutine calc_covariance(b, s, m, n, cb, WX, WY, C)
      integer(IK), intent(in)     :: b(*), s, m, n, cb
      real(RK), intent(in)        :: WX(D, m, n)
      real(RK), intent(inout)     :: WY(D, m), C(cb, n)
      integer(IK)                 :: i, j, ic
!
      ic = 2
      do concurrent(i=1:n)
#ifdef REAL32
        call SGEMM('N', 'T', D, D, m, ONE, WY, D, WX(1, 1, i), D, ZERO, C(ic, i), D)
#else
        call DGEMM('N', 'T', D, D, m, ONE, WY, D, WX(1, 1, i), D, ZERO, C(ic, i), D)
#endif
      end do
!
      do j = 1, s - 1
        ic = ic + DD
        call mol_block_swap(b, j, WY)
        do concurrent(i=1:n)
#ifdef REAL32
          call SGEMM('N', 'T', D, D, m, ONE, WY, D, WX(1, 1, i), D, ZERO, C(ic, i), D)
#else
          call DGEMM('N', 'T', D, D, m, ONE, WY, D, WX(1, 1, i), D, ZERO, C(ic, i), D)
#endif
        end do
        call mol_block_inverse_swap(b, j, WY)
      end do
!
    end subroutine calc_covariance
!
  end subroutine c_matrix_eval
!
!| Calc \( G:=\text{tr}[\mathbf{X}\mathbf{X}^\top]+\text{tr}[\mathbf{Y}\mathbf{Y}^\top] \). <br>
  pure subroutine c_matrix_autocorr(q, C, G)
    integer(IK), intent(in) :: q(*)
    !! c_matrix
    real(RK), intent(in)    :: C(*)
    !! main memory, calculated by c_matrix_eval.
    real(RK), intent(inout) :: G
    !! partial covariance matrix, must be larger than \(d^2\).
    integer(IK)             :: i, nn, nc
    nn = q(nl) * q(nl) * q(cb)
    nc = (1 + q(nl)) * q(cb)
    G = ZERO
    do i = 1, nn, nc
      G = G + C(i)
    end do
  end subroutine c_matrix_autocorr
!
!| Adds \( G_{IJ} \) and \( \mathbf{C}_{IJs} \) specified by index \( i, j, s \) to the arguments. <br>
!  This routine adds directly to G and C(:DD), so they must be initialized. <br>
!  If indices outside the area defined by q(*) is specified, operation results are not guaranteed. <br>
  pure subroutine c_matrix_add(q, i, j, s, C, G, Cp)
    integer(IK), intent(in) :: q(*)
    !! c_matrix
    integer(IK), intent(in) :: i
    !! row index
    integer(IK), intent(in) :: j
    !! collumn index
    integer(IK), intent(in) :: s
    !! symmetry index
    real(RK), intent(in)    :: C(*)
    !! main memory, calculated by c_matrix_eval.
    real(RK), intent(inout) :: G
    !! partial auto variance matrix.
    real(RK), intent(inout) :: Cp(*)
    !! partial covariance matrix, must be larger than \(d^2\).
    integer(IK)             :: k
!
    k = q(cl) * (i - 1) + q(cb) * (j - 1) + 1
    G = G + C(k)
    k = k + 1 + DD * (s - 1)
    call xpy(DD, C(k), Cp)
!
  end subroutine c_matrix_add
!
  pure subroutine zfill(d, x, ld)
    integer(IK), intent(in) :: d
    real(RK), intent(inout) :: x(*)
    integer(IK), intent(in) :: ld
    integer(IK)             :: i, dld
    dld = d * ld
    do concurrent(i=1:dld:ld)
      x(i) = ZERO
    end do
  end subroutine zfill
!
  pure subroutine xpy(N, X, Y)
    integer(IK), intent(in) :: N
    real(RK), intent(in)    :: X(*)
    real(RK), intent(inout) :: Y(*)
    integer(IK)             :: i
    do concurrent(i=1:N)
      Y(i) = X(i) + Y(i)
    end do
  end subroutine xpy
!
  pure function dot(N, X, Y) result(res)
    integer(IK), intent(in) :: N
    real(RK), intent(in)    :: X(*), Y(*)
    real(RK)                :: res
    integer(IK)             :: i
    res = ZERO
    do i = 1, N
      res = res + X(i) * Y(i)
    end do
  end function dot
!
  pure subroutine copy(N, X, Y)
    integer(IK), intent(in) :: N
    real(RK), intent(in)    :: X(*)
    real(RK), intent(inout) :: Y(*)
    integer(IK)             :: i
    do concurrent(i=1:N)
      Y(i) = X(i)
    end do
  end subroutine copy
!
end module mod_c_matrix

