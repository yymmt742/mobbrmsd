      pure function ILAZLR( M, N, A, LDA )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     :: M, N, LDA
      complex(wp),intent(in) :: A( LDA, * )
      integer                :: ILAZLR
      end function ILAZLR

