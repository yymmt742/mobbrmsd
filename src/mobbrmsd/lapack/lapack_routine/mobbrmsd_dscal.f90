!|    mobbrmsd_DSCAL scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!
!  reference DSCAL is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- Reference BLAS level1 routine --
!
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!      jack dongarra, linpack, 3/11/78.
!      modified 3/93 to return if incx .le. 0.
!      modified 12/3/93, array(1) declarations changed to array(*)
!
pure subroutine mobbrmsd_DSCAL(N, DA, DX, INCX)
  integer, intent(in)     :: N
!!         number of elements in input vector(s)
!!
  real(RK), intent(in)    :: DA
!!           On entry, DA specifies the scalar alpha.
!!
  real(RK), intent(inout) :: DX(*)
!!          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!!
  integer, intent(in)     :: INCX
!!         storage spacing between elements of DX
!!
  integer                :: I, M, MP1, NINCX
  intrinsic              :: MOD
!
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

