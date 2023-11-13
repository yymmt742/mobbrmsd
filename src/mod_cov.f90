!| Calculate the covariance matrix.
module mod_cov
  use mod_params, only: IK, RK, ONE => RONE, ZERO => RZERO
  use mod_optarg
  implicit none
  private
  public :: cov, cov_row_major
!
  interface
    include 'dgemm.h'
  end interface
!
  interface cov
    module procedure :: cov_full, cov_part
  end interface cov
!
  interface cov_row_major
    module procedure :: cov_row_major_full, cov_row_major_part
  end interface cov_row_major
!
contains
!
!| Calculate the covariance matrix, X@YT.
  pure subroutine cov_full(d, n, x, y, res, reset)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: n
    !! matrix row dimension.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK), intent(inout) :: res(*)
    !! returns d*d covariance matrix
    logical, intent(in), optional :: reset
    !! reset flag
!
    if (optarg(reset, .true.)) then
      call DGEMM('N', 'T', d, d, n, ONE, x, d, y, d, ZERO, res, d)
    else
      call DGEMM('N', 'T', d, d, n, ONE, x, d, y, d, ONE, res, d)
    end if
!
  end subroutine cov_full
!
!| Calculate the covariance matrix, X@YT, with mask.
  pure subroutine cov_part(d, nlist, x, y, res, reset)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: nlist(:)
    !! matrix row index list.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK), intent(inout) :: res(*)
    !! returns d*d covariance matrix
    logical, intent(in), optional :: reset
    !! reset flag
    integer(IK)             :: i, j, k, p, n
!
    n = SIZE(nlist)
    if (optarg(reset, .true.)) res(:d * d) = ZERO
!
    do i = 1, n
      p = d * (nlist(i) - 1)
      do concurrent(j=1:d, k=1:d)
        block
          integer :: q
          q = (k - 1) * d + j
          res(q) = res(q) + x(p + j) * y(p + k)
        end block
      end do
    end do
!
  end subroutine cov_part
!
!| Calculate the covariance matrix, XT@Y.
  pure subroutine cov_row_major_full(d, n, x, y, res, reset)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: n
    !! matrix row dimension.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK), intent(inout) :: res(*)
    !! returns n*n covariance matrix
    logical, intent(in), optional :: reset
    !! reset flag
!
    if (optarg(reset, .true.)) then
      call DGEMM('T', 'N', n, n, d, ONE, x, d, y, d, ZERO, res, n)
    else
      call DGEMM('T', 'N', n, n, d, ONE, x, d, y, d, ONE, res, n)
    end if
!
  end subroutine cov_row_major_full
!
!| Calculate the covariance matrix, XT@Y, with mask.
  pure subroutine cov_row_major_part(d, nlist, x, y, res, reset)
    integer(IK), intent(in) :: d
    !! matrix collumn dimension.
    integer(IK), intent(in) :: nlist(:)
    !! matrix row index list.
    real(RK), intent(in)    :: x(*)
    !! d*n array
    real(RK), intent(in)    :: y(*)
    !! d*n array
    real(RK), intent(inout) :: res(*)
    !! returns d*d covariance matrix
    logical, intent(in), optional :: reset
    !! reset flag
    integer(IK)             :: i, j, k, n, dd
!
    n = SIZE(nlist)
    dd = d * d
    if (optarg(reset, .true.)) res(:n * n) = ZERO
!
    do concurrent(j=1:n, i=1:n)
      block
        integer :: p, q, r
        p = d * (nlist(i) - 1)
        q = d * (nlist(j) - 1)
        r = i + n * (j - 1)
        do k = 1, d
          res(r) = res(r) + x(k + p) * y(k + q)
        end do
      end block
    end do
!
  end subroutine cov_row_major_part
!
end module mod_cov
