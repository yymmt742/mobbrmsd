      pure subroutine DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
      use LA_CONSTANTS, only: wp=>dp
      CHARACTER(*),intent(in) :: SIDE
      integer,intent(in)      :: LDC, M, N
      real(wp),intent(in)     ::   TAU
      real(wp),intent(in)     :: V( * )
      real(wp),intent(inout)  :: C( LDC, * )
      real(wp),intent(out)    :: WORK( * )
      end subroutine DLARFX
