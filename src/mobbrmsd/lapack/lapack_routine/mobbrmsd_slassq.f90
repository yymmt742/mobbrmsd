!> \brief \b mobbrmsd_SLASSQ updates a sum of squares represented in scaled form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download mobbrmsd_SLASSQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slassq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slassq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slassq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_SLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL               SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       REAL               X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> mobbrmsd_SLASSQ  returns the values  scl  and  smsq  such that
!>
!>    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!>
!> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!> assumed to be non-negative and  scl  returns the value
!>
!>    scl = max( scale, abs( x( i ) ) ).
!>
!> scale and sumsq must be supplied in SCALE and SUMSQ and
!> scl and smsq are overwritten on SCALE and SUMSQ respectively.
!>
!> The routine makes only one pass through the vector x.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements to be used from the vector X.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension (1+(N-1)*INCX)
!>          The vector for which a scaled sum of squares is computed.
!>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of the vector X.
!>          INCX > 0.
!> \endverbatim
!>
!> \param[in,out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          On entry, the value  scale  in the equation above.
!>          On exit, SCALE is overwritten with  scl , the scaling factor
!>          for the sum of squares.
!> \endverbatim
!>
!> \param[in,out] SUMSQ
!> \verbatim
!>          SUMSQ is REAL
!>          On entry, the value  sumsq  in the equation above.
!>          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!>          squares from which  scl  has been factored out.
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
pure subroutine mobbrmsd_SLASSQ(N, X, INCX, SCL, SUMSQ)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
  integer, intent(in) :: INCX, N
  real(RK), intent(inout) :: SCL, SUMSQ
!..
!..Array Arguments..
  real(RK), intent(in)    :: X(*)
!..
!
! =====================================================================
!
!..Local Scalars..
  integer :: IX
  real :: ABSXI
!..
!..Parameters..
! real, parameter :: ZERO = 0.0E+0
!..
! interface
!..external Functions..
!   include 'sisnan.h'
! end interface
!..
!..intrinsic Functions..
  intrinsic :: ABS
!..
!..Executable Statements..
!
  if (N > 0) then
    do IX = 1, 1 + (N - 1) * INCX, INCX
      ABSXI = ABS(X(IX))
      if (ABSXI > ZERO .or. IEEE_IS_NAN(ABSXI)) then
        if (SCL < ABSXI) then
          SUMSQ = 1 + SUMSQ * (SCL / ABSXI)**2
          SCL = ABSXI
        else
          SUMSQ = SUMSQ + (ABSXI / SCL)**2
        end if
      end if
    end do
  end if
  return
!
! end of mobbrmsd_SLASSQ
!
end

