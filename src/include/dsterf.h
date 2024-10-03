      pure subroutine DSTERF( N, D, E, INFO )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     :: N
      integer,intent(out)    :: INFO
      real(wp),intent(inout) :: D( * ), E( * )
      end subroutine DSTERF
