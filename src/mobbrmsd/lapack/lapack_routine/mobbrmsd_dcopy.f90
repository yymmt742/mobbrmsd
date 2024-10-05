!|    mobbrmsd_DCOPY copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to 1.
!
!     Reference IDAMAX is provided by [](http://www.netlib.org/lapack/)
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!      jack dongarra, linpack, 3/11/78.
!      modified 12/3/93, array(1) declarations changed to array(*)
!
pure subroutine mobbrmsd_DCOPY(N, DX, INCX, DY, INCY)
  integer, intent(in)   :: N
!!         number of elements in input vector(s)
  integer, intent(in)   :: INCX
!!         storage spacing between elements of DX
  integer, intent(in)   :: INCY
!!         storage spacing between elements of DY
  real(RK), intent(in)  :: DX(*)
!!          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
  real(RK), intent(out) :: DY(*)
!!          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!
  if (INCX == 1 .and. INCY == 1) then
!
! code for both increments equal to 1
! clean-up loop
!
    block
      integer   :: I, M, MP1
      intrinsic :: MOD
!
      M = MOD(N, 7)
      if (M /= 0) then
        do concurrent(I=1:M)
          DY(I) = DX(I)
        end do
        if (N < 7) return
      end if
!
      MP1 = M + 1
!
      do concurrent(I=MP1:N:7)
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
      end do
!
    end block
!
  elseif ((INCX == 0) .and. (INCY == 0)) then
!
    DY(1) = DX(1)
!
  elseif (INCX == 0) then
!
    block
      integer   :: I, LY, UY
!
      if (INCY < 0) then
        LY = (-N + 1) * INCY + 1
        UY = 1
      else
        LY = 1
        UY = (N - 1) * INCY + 1
      end if
!
      do concurrent(I=LY:UY:INCY)
        DY(I) = DX(1)
      end do
!
    end block
!
  elseif (INCY == 0) then
!
    if (INCX < 0) then
      DY(1) = DX(1)
    else
      DY(1) = DX((N - 1) * INCY + 1)
    end if
!
  else
!
! code for unequal increments or equal increments
! not equal to 1
!
    block
      integer   :: I, IX, IY
!
      IX = 1
      IY = 1
      if (INCX < 0) IX = (-N + 1) * INCX + 1
      if (INCY < 0) IY = (-N + 1) * INCY + 1
!
      do concurrent(I=0:N - 1)
        block
          integer :: IIX, IIY
          IIX = IX + I * INCX
          IIY = IY + I * INCY
          DY(IIY) = DX(IIX)
        end block
      end do
!
    end block
!
  end if
!
! End of mobbrmsd_DCOPY
!
end subroutine mobbrmsd_DCOPY

