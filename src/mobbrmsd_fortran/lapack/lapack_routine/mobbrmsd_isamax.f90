!| mobbrmsd_ISAMAX finds the index of the first element having maximum absolute value.
!
!  reference ISAMAX is provided by [netlib](http://www.netlib.org/lapack/)
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!      jack dongarra, linpack, 3/11/78.
!      modified 3/93 to return if incx .le. 0.
!      modified 12/3/93, array(1) declarations changed to array(*)
!      November 2017
!
pure function mobbrmsd_ISAMAX(N, SX, INCX)
  implicit none
  integer, intent(in)  :: N
!! number of elements in input vector(s)
!!
  real(RK), intent(in) :: SX(*)
!! DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!!
  integer, intent(in)  :: INCX
!! storage spacing between elements of DX
!!
  integer :: mobbrmsd_ISAMAX
!! The index of the first element having maximum absolute value.
!!
  real(RK) :: SMAX
  integer  :: I, IX
  intrinsic :: ABS
!
  mobbrmsd_ISAMAX = 0
  if (N < 1 .or. INCX <= 0) return
  mobbrmsd_ISAMAX = 1
  if (N == 1) return
  if (INCX == 1) then
!
!    code for increment equal to 1
!
    SMAX = ABS(SX(1))
    do I = 2, N
      if (ABS(SX(I)) > SMAX) then
        mobbrmsd_ISAMAX = I
        SMAX = ABS(SX(I))
      end if
    end do
  else
!
!    code for increment not equal to 1
!
    IX = 1
    SMAX = ABS(SX(1))
    IX = IX + INCX
    do I = 2, N
      if (ABS(SX(IX)) > SMAX) then
        mobbrmsd_ISAMAX = I
        SMAX = ABS(SX(IX))
      end if
      IX = IX + INCX
    end do
  end if
  return
end

