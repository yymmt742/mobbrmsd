        pure elemental function IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
        use LA_CONSTANTS, only: sp
        integer,intent(in)      :: IHI, ILO, ISPEC, LWORK, N
        CHARACTER(*),intent(in) :: NAME, OPTS
        integer                 :: IPARMQ
        end function IPARMQ
