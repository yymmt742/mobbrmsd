      pure subroutine DAXPY(N,DA,DX,INCX,DY,INCY)
      use LA_CONSTANTS, only: wp=>dp
      real(wp), intent(in)    :: DA
      integer,  intent(in)    :: INCX,INCY,N
      real(wp), intent(inout) :: DX(*),DY(*)
      end subroutine DAXPY
