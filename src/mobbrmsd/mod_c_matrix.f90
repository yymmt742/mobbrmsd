!| A module for managing c-matrices, tensor of covariance matrices. <br>
!  \( \{\mathbf{C}_{IJs}\}_{IJs} \) is third-order ( \(M\times M\times s \) ) tensor of matrix \( \mathbf{C}_{IJs}\in\mathbb{R}^{d\times d} \),
!  defined by<br>
!  \[ \mathbf{C}_{IJs} = \mathbf{Y}_J \mathbf{Q}_s \mathbf{X}_I^\top \]
!  \( \mathbf{X}_I \) :: \( I \)-th molecule in reference coordinate, \( \mathbf{X} \in\mathbb{R}^{d\times n}\).<br>
!  \( \mathbf{Y}_J \) :: \( J \)-th molecule in target coordinate,    \( \mathbf{Y} \in\mathbb{R}^{d\times n} \).<br>
!  \( \mathbf{Q}_s \) :: Molecular permutation matrix on \( n \). <br>
!  To quickly find the rotation matrix, \( \mathbf{C}_{IJs} \) is stored with autocorrelation \( G_{IJ} \), defined by <br>
!  \[ G_{IJ} = \text{Tr}\left[\mathbf{X}_I\mathbf{X}_I^\top\right] + \text{Tr}\left[\mathbf{Y}_J\mathbf{Y}_J^\top\right] \]
!  @note
!    \( G_{IJ} \) does not change with respect to molecular symmetry permutation index, \( s \). <br>
!    Therefore, data blocks are stored in three-dimensional array C(bs,M,M)
!    with the leading dimension with \( \text{bs}=1 + Sd^2 \). <br>
!    Data blocks are defined by \( \left[ G_{IJ}, \mathbf{C}_{IJ1}, \mathbf{C}_{IJ2}, \dots, \mathbf{C}_{IJS} \right] \) <br>
!  @endnote
module mod_c_matrix
  use mod_dimspec_functions, only: D, DD, compute_cov, covcopy
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
  public :: c_matrix_swap_indices
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
!| C matrix manager.
!   This derived type is mainly used for passing during initialization.
  type c_matrix
    integer(IK) :: q(header_size)
    !! header
    integer(IK), allocatable :: s(:)
    !! state
  contains
    final :: c_matrix_destroy
  end type c_matrix
!
!| Constructer
  interface c_matrix
    module procedure c_matrix_new
  end interface c_matrix
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
    res%q(nw) = MAX( &
              &  mol_block_total_size(b) + mol_block_each_size(b), &
              &  2 * res%q(nl) &
              & )
    allocate (res%s(res%q(nl)))
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
    res = MERGE(q(nw), 0, q(nl) > 0)
  end function c_matrix_worksize
!
!| Evaluation the c-matrix.<br>
!  g-matrix is also calculated at the same time.<br>
!  At the end of the calculation, save the c-matrix C(cb,M,M) to C(\*). <br>
!  workarray W(\*) must be larger than c_matrix_worksize(q). <br>
  pure subroutine c_matrix_eval(q, b, s, X, Y, CX, CY, C, W, sort_by_g)
    integer(IK), intent(in)       :: q(*)
    !! c_matrix header.
    integer(IK), intent(in)       :: b(*)
    !! mol block header.
    integer(IK), intent(inout)    :: s(*)
    !! swap_indice
    real(RK), intent(in)          :: X(*)
    !! reference coordinate,
    !  \(\mathbf{X}\). Arrays must be stored as X(D, n_apm, n_mol).
    real(RK), intent(in)          :: Y(*)
    !! target coordinate,
    ! \(\mathbf{Y}\). Arrays must be stored as Y(D, n_apm, n_mol).
    real(RK), intent(in)          :: CX(*)
    !! centroid of \(\mathbf{X}\).
    real(RK), intent(in)          :: CY(*)
    !! centroid of \(\mathbf{Y}\).
    real(RK), intent(inout)       :: C(*)
    !! main memory
    real(RK), intent(inout)       :: W(*)
    !! work memory
    logical, intent(in), optional :: sort_by_g
    !! if true, row is sorted respect to G of reference coordinate.
    integer(IK), parameter     :: gx = 1
    integer(IK), parameter     :: wx = 1
    integer(IK)                :: n_sym, n_apm, gy, wy, i
    associate (n_blk => q(cb), n_mol => q(nl))
      n_apm = mol_block_napm(b)
      do concurrent(i=1:n_mol)
        s(i) = i
      end do
      if (n_apm < 1) then
        call zfill(q(cl) * q(nl), C, 1)
        return
      end if
      n_sym = mol_block_nsym(b)
      gy = gx + n_mol
      wy = wx + mol_block_total_size(b)
      call eval_g(D, n_apm, n_mol, X, CX, W(gx))
      call eval_g(D, n_apm, n_mol, Y, CY, W(gy))
      if (PRESENT(sort_by_g)) then
        if (sort_by_g .and. n_mol > 1) call qs(n_mol, s, W(gx))
      else
        if (n_mol > 1) call qs(n_mol, s, W(gx))
      end if
      call set_g(n_mol, n_blk, W(gx), W(gy), C)
      call eval_c_matrix(b, s, n_sym, n_apm, n_mol, n_blk, X, Y, CX, CY, C, W(wx), W(wy))
    end associate
  end subroutine c_matrix_eval
!
! trace of self correlation matrix
  pure subroutine eval_g(D, n_apm, n_mol, X, CX, G)
    integer(IK), intent(in)        :: D, n_apm, n_mol
    real(RK), intent(in)           :: X(D, n_apm, n_mol), CX(D)
    real(RK), intent(inout)        :: G(n_mol)
    integer(IK)                    :: i, j, k
    do concurrent(k=1:n_mol)
      G(k) = ZERO
      do j = 1, n_apm
        do i = 1, D
          G(k) = G(k) + (X(i, j, k) - CX(i))**2
        end do
      end do
    end do
  end subroutine eval_g
!
! trace of self correlation matrix
  pure subroutine set_g(n_mol, n_blk, GX, GY, C)
    integer(IK), intent(in)        :: n_mol, n_blk
    real(RK), intent(in)           :: GX(*), GY(*)
    real(RK), intent(inout)        :: C(n_blk, n_mol, n_mol)
    integer(IK)                    :: i, j
    do concurrent(i=1:n_mol, j=1:n_mol)
      C(1, i, j) = GX(i) + GY(j)
    end do
  end subroutine set_g
!
! get correlation matrix \(\mathbf{C} = \mathbf{Y}^{\top}\mathbf{X}\)
! and optimal rotation \(\mathbf{R}^{\top}\)
  pure subroutine eval_c_matrix(b, s, n_sym, n_apm, n_mol, n_blk, X, Y, CX, CY, C, WX, WY)
    integer(IK), intent(in)     :: b(*), s(*), n_sym, n_apm, n_mol, n_blk
    real(RK), intent(in)        :: X(D, n_apm, n_mol), Y(D, n_apm, n_mol)
    real(RK), intent(in)        :: CX(D), CY(D)
    real(RK), intent(inout)     :: C(n_blk, n_mol, n_mol), WX(D, n_apm, n_mol), WY(D, n_apm)
    integer(IK)                 :: i
    do concurrent(i=1:n_mol)
      call covcopy(D, n_apm, X(1, 1, s(i)), CX, WX(1, 1, i))
    end do
    do i = 1, n_mol
      call covcopy(D, n_apm, Y(1, 1, i), CY, WY)
      call calc_covariance(b, n_sym, n_apm, n_mol, n_blk, WX, WY, C(1, 1, i))
    end do
  end subroutine eval_c_matrix
!
  pure subroutine calc_covariance(b, n_sym, n_apm, n_mol, n_blk, WX, WY, C)
    integer(IK), intent(in)     :: b(*), n_sym, n_apm, n_mol, n_blk
    real(RK), intent(in)        :: WX(D, n_apm, n_mol)
    real(RK), intent(inout)     :: WY(D, n_apm), C(n_blk, n_mol)
    integer(IK)                 :: i, j, ic
    ic = 2
    do concurrent(i=1:n_mol)
      call compute_cov(D, n_apm, WX(1, 1, i), WY, C(ic, i))
    end do
    do j = 1, n_sym - 1
      ic = ic + DD
      call mol_block_swap(b, j, WY)
      do concurrent(i=1:n_mol)
        call compute_cov(D, n_apm, WX(1, 1, i), WY, C(ic, i))
      end do
      call mol_block_inverse_swap(b, j, WY)
    end do
  end subroutine calc_covariance
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
!| returns swap(z)
  pure subroutine c_matrix_swap_indices(q, s, z, res)
    integer(IK), intent(in)    :: q(*)
    !! header
    integer(IK), intent(in)    :: s(*)
    !! state
    integer(IK), intent(in)    :: z(*)
    !! permutation
    integer(IK), intent(inout) :: res(*)
    !! swap indice
    integer(IK)                :: i
    do concurrent(i=1:q(nl))
      res(i) = s(z(i))
    end do
  end subroutine c_matrix_swap_indices
!
  pure elemental subroutine c_matrix_destroy(this)
    type(c_matrix), intent(inout) :: this
    if (ALLOCATED(this%s)) deallocate (this%s)
  end subroutine c_matrix_destroy
!
! ---
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
! pure subroutine copy(N, X, Y)
!   integer(IK), intent(in) :: N
!   real(RK), intent(in)    :: X(*)
!   real(RK), intent(inout) :: Y(*)
!   integer(IK)             :: i
!   do concurrent(i=1:N)
!     Y(i) = X(i)
!   end do
! end subroutine copy
!
!| Quick sort
  pure recursive subroutine qs(n, s, g)
    integer, intent(in)     :: n
    integer, intent(inout)  :: s(*)
    real(RK), intent(inout) :: g(*)
    real(RK)                :: h
    integer                 :: i, j, p, t
    j = n; i = 1; p = j / 2
    do
      do while (g(i) > g(p)); i = i + 1; end do
      do while (g(j) < g(p)); j = j - 1; end do
      if (i >= j) exit
      h = g(i); g(i) = g(j); g(j) = h
      t = s(i); s(i) = s(j); s(j) = t
      if (i == p) then; p = j
      elseif (j == p) then; p = i
      end if
      i = i + 1; j = j - 1
    end do
    if (2 < i) call qs(i - 1, s(1), g(1))
    if (j < n - 1) call qs(n - j, s(j + 1), g(j + 1))
  end subroutine qs
!
end module mod_c_matrix

