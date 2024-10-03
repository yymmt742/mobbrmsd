      pure function DLANST( NORM, N, D, E )
      use LA_CONSTANTS, only: wp=>dp
      character(*),intent(in) :: NORM
      integer,intent(in)      :: N
      real(wp),intent(in)     :: D( * ), E( * )
      real(wp)                :: DLANST
      end function DLANST
