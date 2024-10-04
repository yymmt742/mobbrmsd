! > \brief\b mobbrmsd_ILASLR scans a matrix for its last non - zero row.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
! http://www.netlib.org/lapack/explore-html/
!
! > \htmlonly
! > Download mobbrmsd_ILASLR + dependencies
! >  < a href = "http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaslr.f" >
! > [TGZ] < /a >
! >  < a href = "http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaslr.f" >
! > [ZIP] < /a >
! >  < a href = "http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaslr.f" >
! > [TXT] < /a >
! > \endhtmlonly
!
!  Definition:
!  ===========
!
!   INTEGER FUNCTION mobbrmsd_ILASLR( M, N, A, LDA )
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
! > mobbrmsd_ILASLR scans A for its last non - zero row.
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
! > \date December 2016
!
! > \ingroup realOTHERauxiliary
!
!  =====================================================================
pure function mobbrmsd_ILASLR(M, N, A, LDA)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! December 2016
!
! .. Scalar Arguments ..
  integer, intent(in)  :: M, N, LDA
! ..
! .. Array Arguments ..
  real(RK), intent(in) :: A(LDA, *)
  integer :: mobbrmsd_ILASLR
! ..
!  =====================================================================
! ..
! .. Local Scalars ..
  integer :: I, J
! ..
! .. Executable Statements ..
!
! Quick test for the common case where one corner is non-zero.
  if (M == 0) then
    mobbrmsd_ILASLR = M
  elseif (A(M, 1) /= ZERO .or. A(M, N) /= ZERO) then
    mobbrmsd_ILASLR = M
  else
! Scan up each column tracking the last zero row seen.
    mobbrmsd_ILASLR = 0
    do J = 1, N
      I = M
      do while ((A(MAX(I, 1), J) == ZERO) .and. (I >= 1))
        I = I - 1
      end do
      mobbrmsd_ILASLR = MAX(mobbrmsd_ILASLR, I)
    end do
  end if
  return
end
