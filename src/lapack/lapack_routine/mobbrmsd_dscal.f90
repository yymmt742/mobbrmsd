!> \brief \b mobbrmsd_DSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE mobbrmsd_DSCAL(N,DA,DX,INCX)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DA
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
!>    mobbrmsd_DSCAL scales a vector by a constant.
!>    uses unrolled loops for increment equal to 1.
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
!> \param[in] DA
!> \verbatim
!>          DA is DOUBLE PRECISION
!>           On entry, DA specifies the scalar alpha.
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
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
pure subroutine mobbrmsd_DSCAL(N, DA, DX, INCX)
! use LA_CONSTANTS, only: RK => dp
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
  integer, intent(in)     :: INCX, N
  real(RK), intent(in)    :: DA
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: DX(*)
!     ..
!  =====================================================================
!     .. Local Scalars ..
  integer                :: I, M, MP1, NINCX
!     ..
!     .. Intrinsic Functions ..
  intrinsic              :: MOD
!     ..
  if (N <= 0 .or. INCX <= 0) return
  if (INCX == 1) then
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
    M = MOD(N, 5)
    if (M /= 0) then
      do I = 1, M
        DX(I) = DA * DX(I)
      end do
      if (N < 5) return
    end if
    MP1 = M + 1
    do I = MP1, N, 5
      DX(I) = DA * DX(I)
      DX(I + 1) = DA * DX(I + 1)
      DX(I + 2) = DA * DX(I + 2)
      DX(I + 3) = DA * DX(I + 3)
      DX(I + 4) = DA * DX(I + 4)
    end do
  else
!
!        code for increment not equal to 1
!
    NINCX = N * INCX
    do I = 1, NINCX, INCX
      DX(I) = DA * DX(I)
    end do
  end if
  return
!
!     End of mobbrmsd_DSCAL
!
end subroutine mobbrmsd_DSCAL
