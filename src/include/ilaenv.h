        pure function ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
        character( * ),intent(in) :: NAME, OPTS
        integer,intent(in)        :: ISPEC, N1, N2, N3, N4
        integer                   :: ILAENV
        end function ILAENV
