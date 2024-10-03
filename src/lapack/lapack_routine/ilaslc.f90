! > \brief\b ILASLC scans a matrix for its last non - zero column.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
! http://www.netlib.org/lapack/explore-html/
!
! > \htmlonly
! > Download ILASLC + dependencies
! >  < a href = "http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaslc.f" >
! > [TGZ] < /a >
! >  < a href = "http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaslc.f" >
! > [ZIP] < /a >
! >  < a href = "http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaslc.f" >
! > [TXT] < /a >
! > \endhtmlonly
!
!  Definition:
!  ===========
!
!   INTEGER FUNCTION ILASLC( M, N, A, LDA )
!
!   .. Scalar Arguments ..
!   INTEGER            M, N, LDA
!   ..
!   .. Array Arguments ..
!   REAL               A( LDA, ! )
!   ..
!
!
! > \par Purpose:
!  =============
! >
! > \verbatim
! >
! > ILASLC scans A for its last non - zero column.
! > \endverbatim
!
!  Arguments:
!  ==========
!
! > \param[in] M
! > \verbatim
! > M is integer
! > The number of rows of the matrix A.
! > \endverbatim
! >
! > \param[in] N
! > \verbatim
! > N is integer
! > The number of columns of the matrix A.
! > \endverbatim
! >
! > \param[in] A
! > \verbatim
! > A is real array, dimension(LDA, N)
! > The m by n matrix A.
! > \endverbatim
! >
! > \param[in] LDA
! > \verbatim
! > LDA is integer
! > The leading dimension of the array A.LDA >= MAX(1, M) .
! > \endverbatim
!
!  Authors:
!  ========
!
! > \author Univ.of Tennessee
! > \author Univ.of California Berkeley
! > \author Univ.of Colorado Denver
! > \author NAG Ltd.
!
! > \date June 2017
!
! > \ingroup realOTHERauxiliary
!
!  =====================================================================
pure function ILASLC(M, N, A, LDA)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! June 2017
!
! .. Scalar Arguments ..
  integer, intent(in)  :: M, N, LDA
! ..
! .. Array Arguments ..
  real(RK), intent(in) :: A(LDA, *)
  integer              :: ILASLC
! ..
!  =====================================================================
!
! .. Local Scalars ..
  integer :: I
! ..
! .. Executable Statements ..
!
! Quick test for the common case where one corner is non-zero.
  if (N == 0) then
    ILASLC = N
  else if (A(1, N) /= ZERO .or. A(M, N) /= ZERO) then
    ILASLC = N
  else
! Now scan each column from the end, returning with the first non-zero.
    do ILASLC = N, 1, -1
      do I = 1, M
        if (A(I, ILASLC) /= ZERO) return
      end do
    end do
  end if
  return
end
