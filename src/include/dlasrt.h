        pure subroutine DLASRT( ID, N, D, INFO )
        use LA_CONSTANTS, only: wp=>dp
        character(*),intent(in) :: ID
        integer,intent(in)      :: N
        integer,intent(out)     :: INFO
        real(wp),intent(inout)  :: D( * )
        end subroutine DLASRT
