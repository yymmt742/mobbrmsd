!> \brief \b mobbrmsd_IDAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION mobbrmsd_IDAMAX(N,DX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    mobbrmsd_IDAMAX finds the index of the first element having maximum absolute value.
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
!> \param[in] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
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
!> \ingroup aux_blas
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
pure function mobbrmsd_IDAMAX(N, DX, INCX)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)  :: INCX, N
!     ..
!     .. Array Arguments ..
  real(RK), intent(in) :: DX(*)
!     ..
!
  integer :: mobbrmsd_IDAMAX
!
!  =====================================================================
!
!     .. Local Scalars ..
  real(RK) :: DMAX
  integer :: I, IX
!     ..
!     .. Intrinsic Functions ..
  intrinsic :: ABS
!     ..
  mobbrmsd_IDAMAX = 0
  if (N < 1 .or. INCX <= 0) return
  mobbrmsd_IDAMAX = 1
  if (N == 1) return
  if (INCX == 1) then
!
!        code for increment equal to 1
!
    DMAX = ABS(DX(1))
    do I = 2, N
      if (ABS(DX(I)) > DMAX) then
        mobbrmsd_IDAMAX = I
        DMAX = ABS(DX(I))
      end if
    end do
  else
!
!        code for increment not equal to 1
!
    IX = 1
    DMAX = ABS(DX(1))
    IX = IX + INCX
    do I = 2, N
      if (ABS(DX(IX)) > DMAX) then
        mobbrmsd_IDAMAX = I
        DMAX = ABS(DX(IX))
      end if
      IX = IX + INCX
    end do
  end if
  return
!
!     End of mobbrmsd_IDAMAX
!
end function mobbrmsd_IDAMAX

