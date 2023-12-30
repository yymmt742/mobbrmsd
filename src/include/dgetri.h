pure subroutine DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
  use mod_params, only : wp=>R8
  integer,intent(in)     :: LDA, LWORK, N
  integer,intent(out)    :: INFO
  integer,intent(in)     :: IPIV( * )
  real(wp),intent(inout) :: A( LDA, * )
  real(wp),intent(out)   :: WORK( * )
end subroutine DGETRI
