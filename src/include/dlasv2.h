pure subroutine DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
use LA_CONSTANTS, only: DP
real(DP),intent(in)  :: F, G, H
real(DP),intent(out) :: SSMAX, SSMIN, SNL, SNR, CSL, CSR
end subroutine DLASV2
