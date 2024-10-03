      pure function ILADLR( M, N, A, LDA )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)  :: M, N, LDA
      real(wp),intent(in) :: A( LDA, * )
      integer             :: ILADLR
      end function ILADLR

