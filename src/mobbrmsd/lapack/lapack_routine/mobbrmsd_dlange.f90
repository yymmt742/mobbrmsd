!> \brief \b mobbrmsd_DLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_DLANGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlange.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlange.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlange.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION mobbrmsd_DLANGE( NORM, M, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_DLANGE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> real matrix A.
!> \endverbatim
!>
!> \return mobbrmsd_DLANGE
!> \verbatim
!>
!>    mobbrmsd_DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in mobbrmsd_DLANGE as described
!>          above.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.  When M = 0,
!>          mobbrmsd_DLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.  When N = 0,
!>          mobbrmsd_DLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(M,1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!>          referenced.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup doubleGEauxiliary
!
!  =====================================================================
pure subroutine mobbrmsd_DLANGE(NORM, M, N, A, LDA, RES, WORK)
! use LA_CONSTANTS, only: RK => dp
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  character(*), intent(in) :: NORM
  integer, intent(in)      :: LDA, M, N
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout)  :: A(LDA, *)
  real(RK), intent(out)    :: RES, WORK(*)
!     ..
!
! =====================================================================
!     ..
!     .. Local Scalars ..
  integer                  :: I, J
  double precision         :: SCALE, SUM, TEMP
!     ..
!     .. Parameters ..
! real(RK), parameter      :: ONE = 1.0_RK
! real(RK), parameter      :: ZERO = 0.0_RK
! interface
!     .. External Subroutines ..
!   include 'dlassq.h'
!     .. External Functions ..
!   include 'lsame.h'
!   include 'disnan.h'
! end interface
!     ..
!     .. Intrinsic Functions ..
  intrinsic                 :: ABS, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
  if (MIN(M, N) == 0) then
    RES = ZERO
  else if (mobbrmsd_LSAME(NORM, 'M')) then
!
!        Find max(abs(A(i,j))).
!
    RES = ZERO
    do J = 1, N
      do I = 1, M
        TEMP = ABS(A(I, J))
        if (RES < TEMP .or. IEEE_IS_NAN(TEMP)) RES = TEMP
      end do
    end do
  else if ((mobbrmsd_LSAME(NORM, 'O')) .or. (NORM == '1')) then
!
!         Find norm1(A).
!
    RES = ZERO
    do J = 1, N
      SUM = ZERO
      do I = 1, M
        SUM = SUM + ABS(A(I, J))
      end do
      if (RES < SUM .or. IEEE_IS_NAN(SUM)) RES = SUM
    end do
  else if (mobbrmsd_LSAME(NORM, 'I')) then
!
!         Find normI(A).
!
    do I = 1, M
      WORK(I) = ZERO
    end do
    do J = 1, N
      do I = 1, M
        WORK(I) = WORK(I) + ABS(A(I, J))
      end do
    end do
    RES = ZERO
    do I = 1, M
      TEMP = WORK(I)
      if (RES < TEMP .or. IEEE_IS_NAN(TEMP)) RES = TEMP
    end do
  else if ((mobbrmsd_LSAME(NORM, 'F')) .or. (mobbrmsd_LSAME(NORM, 'E'))) then
!
!  Find normF(A).
!
    SCALE = ZERO
    SUM = ONE
    do J = 1, N
      call mobbrmsd_DLASSQ(M, A(1, J), 1, SCALE, SUM)
    end do
    RES = SCALE * SQRT(SUM)
  end if
!
  return
!
!     End of mobbrmsd_DLANGE
!
end subroutine mobbrmsd_DLANGE

