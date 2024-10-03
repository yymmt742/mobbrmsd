!> \brief \b mobbrmsd_SNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION mobbrmsd_SNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_SNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    mobbrmsd_SNRM2 := sqrt( x'*x ).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
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
!> \date November 2017
!
!> \ingroup single_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  -- This version written on 25-October-1982.
!>     Modified on 14-October-1993 to inline the call to mobbrmsd_SLASSQ.
!>     Sven Hammarling, Nag Ltd.
!> \endverbatim
!>
!  =====================================================================
pure function mobbrmsd_SNRM2(N, X, INCX)
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
  integer, intent(in) :: INCX, N
!..
!..Array Arguments..
  real(RK), intent(in) :: X(*)
  real(RK)             :: mobbrmsd_SNRM2
!..
!
!  =====================================================================
!
!..Parameters..
! real(RK), parameter :: ONE = 1.0E+0, ZERO = 0.0E+0
!..
!..Local Scalars..
  real(RK) :: ABSXI, NORM, SCL, SSQ
  integer :: IX
!..
!..intrinsic Functions..
  intrinsic :: ABS, SQRT
!..
  if (N < 1 .or. INCX < 1) then
    NORM = ZERO
  else if (N == 1) then
    NORM = ABS(X(1))
  else
    SCL = ZERO
    SSQ = ONE
! The following loop is equivalent to this call to the LAPACK
! auxiliary routine:
! call mobbrmsd_SLASSQ(N, X, INCX, SCL, SSQ)
!
    do IX = 1, 1 + (N - 1) * INCX, INCX
      if (X(IX) /= ZERO) then
        ABSXI = ABS(X(IX))
        if (SCL < ABSXI) then
          SSQ = ONE + SSQ * (SCL / ABSXI)**2
          SCL = ABSXI
        else
          SSQ = SSQ + (ABSXI / SCL)**2
        end if
      end if
    end do
    NORM = SCL * SQRT(SSQ)
  end if
!
  mobbrmsd_SNRM2 = NORM
  return
!
! end of mobbrmsd_SNRM2.
!
end
