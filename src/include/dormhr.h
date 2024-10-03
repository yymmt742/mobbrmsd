      pure subroutine DORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C, &
     &                   LDC, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: SIDE, TRANS
      integer,intent(in)      :: IHI, ILO, LDA, LDC, LWORK, M, N
      integer,intent(out)     :: INFO
      real(wp),intent(in)     :: TAU( * )
      real(wp),intent(inout)  :: A( LDA, * ), C( LDC, * )
      real(wp),intent(out)    :: WORK( * )
      end subroutine DORMHR
