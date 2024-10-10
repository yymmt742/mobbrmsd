!|    mobbrmsd_SSCAL scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!
!  Reference SSCAL is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!      jack dongarra, linpack, 3/11/78.
!      modified 3/93 to return if incx .le. 0.
!      modified 12/3/93, array(1) declarations changed to array(*)
!
pure subroutine mobbrmsd_SSCAL(N, SA, SX, INCX)
  implicit none
  integer, intent(in)  :: N
!!         number of elements in input vector(s)
!!
  integer, intent(in)  :: INCX
!!         storage spacing between elements of SX
!!
  real(RK), intent(in) :: SA
!!           On entry, SA specifies the scalar alpha.
!!
  real(RK), intent(inout) :: SX(*)
!!          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
  integer :: I, M, MP1, NINCX
  intrinsic :: MOD
!
  if (N <= 0 .or. INCX <= 0) return
  if (INCX == 1) then
!
! code for increment equal to 1
!
!
! clean - up loop
!
    M = MOD(N, 5)
    if (M /= 0) then
      do I = 1, M
        SX(I) = SA * SX(I)
      end do
      if (N < 5) return
    end if
    MP1 = M + 1
    do I = MP1, N, 5
      SX(I) = SA * SX(I)
      SX(I + 1) = SA * SX(I + 1)
      SX(I + 2) = SA * SX(I + 2)
      SX(I + 3) = SA * SX(I + 3)
      SX(I + 4) = SA * SX(I + 4)
    end do
  else
!
! code for increment not equal to 1
!
    NINCX = N * INCX
    do I = 1, NINCX, INCX
      SX(I) = SA * SX(I)
    end do
  end if
  return
end

