        pure subroutine DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, &
       &                        DN1, DN2, TAU, TTYPE, G )
        use LA_CONSTANTS, only: wp=>dp
        integer,intent(in)     :: I0, N0, PP, N0IN
        integer,intent(out)    :: TTYPE
        real(wp),intent(in)    :: DMIN, DMIN1, DMIN2, DN, DN1, DN2
        real(wp),intent(inout) :: G
        real(wp),intent(out)   :: TAU
        real(wp),intent(in)    :: Z( * )
        end subroutine DLASQ4
