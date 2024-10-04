!> \brief \b mobbrmsd_DSWAP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DSWAP(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    mobbrmsd_DSWAP interchanges two vectors.
!>    uses unrolled loops for increments equal to 1.
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
!> \param[in,out] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[in,out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
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
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
pure subroutine mobbrmsd_DSWAP(N, DX, INCX, DY, INCY)
! use LA_CONSTANTS, only: RK => dp
  implicit none
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)     :: INCX, INCY, N
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: DX(*), DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  integer                :: I, IX, IY, M, MP1
!     ..
!     .. Intrinsic Functions ..
  intrinsic              :: MOD
!     ..
  if (N <= 0) return
  if (INCX == 1 .and. INCY == 1) then
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
    M = MOD(N, 3)
    if (M /= 0) then
      do concurrent(I=1:M)
        block
          real(RK) :: DTEMP
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
        end block
      end do
      if (N < 3) return
    end if
    MP1 = M + 1
    do concurrent(I=MP1:N:3)
      block
        real(RK) :: DTEMP
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
      end block
    end do
  else
!
!       code for unequal increments or equal increments not equal
!         to 1
!
    IX = 1
    IY = 1
    if (INCX < 0) IX = (-N + 1) * INCX + 1
    if (INCY < 0) IY = (-N + 1) * INCY + 1
!
    do concurrent(I=0:N - 1)
      block
        real(RK) :: DTEMP
        integer  :: IIX, IIY
        IIX = IX + I * INCX
        IIY = IY + I * INCY
        DTEMP = DX(IIX)
        DX(IIX) = DY(IIY)
        DY(IIX) = DTEMP
      end block
    end do
  end if
  return
!
!     End of mobbrmsd_DSWAP
!
end subroutine mobbrmsd_DSWAP

