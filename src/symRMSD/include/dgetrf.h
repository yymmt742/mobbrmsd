pure subroutine DGETRF( M, N, A, LDA, IPIV, INFO )
  use mod_params, only :  wp=>R8
  integer,intent(in)      :: LDA, M, N
  integer,intent(out)     :: INFO
  integer,intent(out)     :: IPIV( * )
  real(wp),intent(inout)  :: A( LDA, * )
end subroutine DGETRF
