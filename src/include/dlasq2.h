      pure subroutine DLASQ2( N, Z, INFO )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)     ::  N
      integer,intent(out)    ::  INFO
      real(wp),intent(inout) :: Z( * )
      end subroutine DLASQ2
