        pure subroutine DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, &
       &                        DNM1, DNM2 )
        use LA_CONSTANTS, only: wp=>dp
        integer,intent(in)     :: I0, N0, PP
        real(wp),intent(out)   :: DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
        real(wp),intent(inout) :: Z( * )
        end  subroutine DLASQ6
