      pure subroutine DLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )
      use LA_CONSTANTS, only: wp=>dp
      real(wp),intent(in)     :: SI1, SI2, SR1, SR2
      integer,intent(in)      :: LDH, N
      real(wp),intent(in)     :: H( LDH, * )
      real(wp),intent(out)    :: V( * )
      end subroutine DLAQR1
