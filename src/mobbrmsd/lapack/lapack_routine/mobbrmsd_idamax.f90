!| mobbrmsd_IDAMAX finds the index of the first element having maximum absolute value.
!
!  reference IDAMAX is provided by [netlib](http://www.netlib.org/lapack/)
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!      jack dongarra, linpack, 3/11/78.
!      modified 3/93 to return if incx .le. 0.
!      modified 12/3/93, array(1) declarations changed to array(*)
!
pure function mobbrmsd_IDAMAX(N, DX, INCX)
  implicit none
  integer, intent(in)  :: N
!! number of elements in input vector(s)
!!
  real(RK), intent(in) :: DX(*)
!! DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!!
  integer, intent(in)  :: INCX
!! storage spacing between elements of DX
!!
  integer :: mobbrmsd_IDAMAX
!! The index of the first element having maximum absolute value.
!!
  real(RK) :: DMAX
  integer :: I, IX
  intrinsic :: ABS
!
  mobbrmsd_IDAMAX = 0
  if (N < 1 .or. INCX <= 0) return
  mobbrmsd_IDAMAX = 1
  if (N == 1) return
  if (INCX == 1) then
!
! code for increment equal to 1
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
! code for increment not equal to 1
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
! End of mobbrmsd_IDAMAX
!
end function mobbrmsd_IDAMAX

