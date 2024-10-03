        pure subroutine DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,&
       &                   ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,     &
       &                   DN2, G, TAU )
        use LA_CONSTANTS, only: wp=>dp
        integer,intent(in)     :: I0
        integer,intent(inout)  :: N0, PP, NFAIL, ITER, NDIV,TTYPE
        real(wp),intent(out)   :: DMIN,  SIGMA
        real(wp),intent(inout) :: DESIG, QMAX, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
        logical,intent(in)     :: IEEE
        real(wp),intent(inout) :: Z( * )
        end subroutine DLASQ3
