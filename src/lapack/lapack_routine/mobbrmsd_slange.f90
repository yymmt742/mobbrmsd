!> \brief \b mobbrmsd_SLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute valie of any element of a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_SLANGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slange.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slange.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slange.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION mobbrmsd_SLANGE( NORM, M, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_SLANGE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> real matrix A.
!> \endverbatim
!>
!> \return mobbrmsd_SLANGE
!> \verbatim
!>
!>    mobbrmsd_SLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!>          Specifies the value to be returned in mobbrmsd_SLANGE as described
!>          above.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.  When M = 0,
!>          mobbrmsd_SLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.  When N = 0,
!>          mobbrmsd_SLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
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
!>          WORK is REAL array, dimension (MAX(1,LWORK)),
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
!> \date December 2016
!
!> \ingroup realGEauxiliary
!
!  =====================================================================
pure subroutine mobbrmsd_SLANGE(NORM, M, N, A, LDA, RES, WORK)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
  implicit none
!..Scalar Arguments..
  character, intent(in) :: NORM
  integer, intent(in)   :: LDA, M, N
!..
!..Array Arguments..
  real(RK), intent(inout) :: A(LDA, *)
  real(RK), intent(out)   :: WORK(*), RES
!..
!
! =====================================================================
!..
!..Local Scalars..
  integer :: I, J
  real(RK) :: SUM, val, TEMP
!..
!..Local Arrays..
  real(RK) :: SSQ(2), COLSSQ(2)
!
!..Parameters..
! real(RK), parameter :: ONE = 1.0E+0, ZERO = 0.0E+0
!..
! interface
! .. External Functions ..
!   include 'lsame.h'
!   include 'sisnan.h'
! .. External Subroutines ..
!   include 'slassq.h'
!   include 'scombssq.h'
! end interface
!..
!..intrinsic Functions..
  intrinsic :: ABS, MIN, SQRT
!..
!..Executable Statements..
!
  if (MIN(M, N) == 0) then
    val = ZERO
  else if (mobbrmsd_LSAME(NORM, 'M')) then
    !
    !Find MAX(ABS(A(i, j))) .
    !
    val = ZERO
    do J = 1, N
      do I = 1, M
        TEMP = ABS(A(I, J))
        if (val < TEMP .or. mobbrmsd_SISNAN(TEMP)) val = TEMP
      end do
    end do
  else if ((mobbrmsd_LSAME(NORM, 'O')) .or. (NORM == '1')) then
    !
    !Find norm1(A) .
    !
    val = ZERO
    do J = 1, N
      SUM = ZERO
      do I = 1, M
        SUM = SUM + ABS(A(I, J))
      end do
      if (val < SUM .or. mobbrmsd_SISNAN(SUM)) val = SUM
    end do
  else if (mobbrmsd_LSAME(NORM, 'I')) then
    !
    ! Find normI(A) .
    !
    do I = 1, M
      WORK(I) = ZERO
    end do
    do J = 1, N
      do I = 1, M
        WORK(I) = WORK(I) + ABS(A(I, J))
      end do
    end do
    val = ZERO
    do I = 1, M
      TEMP = WORK(I)
      if (val < TEMP .or. mobbrmsd_SISNAN(TEMP)) val = TEMP
    end do
  else if ((mobbrmsd_LSAME(NORM, 'F')) .or. (mobbrmsd_LSAME(NORM, 'E'))) then
    !
    !Find normF(A) .
    !SSQ(1) is scale
    !SSQ(2) is sum - of - squares
    !For better accuracy, sum each column separately.
    !
    SSQ(1) = ZERO
    SSQ(2) = ONE
    do J = 1, N
      COLSSQ(1) = ZERO
      COLSSQ(2) = ONE
      call mobbrmsd_SLASSQ(M, A(1, J), 1, COLSSQ(1), COLSSQ(2))
      call mobbrmsd_SCOMBSSQ(SSQ, COLSSQ)
    end do
    val = SSQ(1) * SQRT(SSQ(2))
  end if
  !
  RES = val
  return
  !
  !end of mobbrmsd_SLANGE
  !
end
