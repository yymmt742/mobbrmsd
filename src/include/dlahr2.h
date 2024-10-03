      pure subroutine DLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     :: K, LDA, LDT, LDY, N, NB
      real(wp),intent(inout) :: A( LDA, * )
      real(wp),intent(out)   :: T( LDT, NB ), TAU( NB ), Y( LDY, NB )
      end subroutine DLAHR2
