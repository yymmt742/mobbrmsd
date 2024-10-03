      pure subroutine DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: COMPZ
      integer,intent(in)      :: LDZ, N
      integer,intent(out)     :: INFO
      real(wp),intent(inout)  :: D( * ), E( * ), Z( LDZ, * )
      real(wp),intent(out)    :: WORK( * )
      end subroutine DSTEQR
