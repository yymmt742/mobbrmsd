!| mobbrmsd_SLASQ6 computes one dqd transform in ping-pong form.
!
!  mobbrmsd_SLASQ6 computes one dqd (shift equal to zero) transform
!  in ping-pong form, with protection against underflow and overflow.
!
!  Reference SLASQ6 is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK computational routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
pure subroutine mobbrmsd_SLASQ6(I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2)
  implicit none
  integer, intent(in)     :: I0
!!  First index.
!!
  integer, intent(in)     :: N0
!!  Last index.
!!
  real(RK), intent(inout) :: Z(*)
!!  DOUBLE PRECISION array, dimension ( 4*N )
!!  Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
!!  an extra argument.
!!
  integer, intent(in)     :: PP
!!  PP=0 for ping, PP=1 for pong.
!!
  real(RK), intent(out)   :: DMIN
!!  Minimum value of d.
!!
  real(RK), intent(out)   :: DMIN1
!!  Minimum value of d, excluding D( N0 ).
!!
  real(RK), intent(out)   :: DMIN2
!!  Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!!
  real(RK), intent(out)   :: DN
!!  d(N0), the last value of d.
!!
  real(RK), intent(out)   :: DNM1
!!  d(N0-1).
!!
  real(RK), intent(out)   :: DNM2
!!  d(N0-2).
!!
  integer :: J4, J4P2
  real(RK) :: D, EMIN, SAFMIN, TEMP
  intrinsic :: MIN
! interface
!   include 'slamch.h'
! end interface
!
  if ((N0 - I0 - 1) <= 0) return
!
  SAFMIN = mobbrmsd_SLAMCH('Safe minimum')
  J4 = 4 * I0 + PP - 3
  EMIN = Z(J4 + 4)
  D = Z(J4)
  DMIN = D
!
  if (PP == 0) then
    do J4 = 4 * I0, 4 * (N0 - 3), 4
      Z(J4 - 2) = D + Z(J4 - 1)
      if (Z(J4 - 2) == ZERO) then
        Z(J4) = ZERO
        D = Z(J4 + 1)
        DMIN = D
        EMIN = ZERO
      else if (SAFMIN * Z(J4 + 1) < Z(J4 - 2) .and. SAFMIN * Z(J4 - 2) < Z(J4 + 1)) then
        TEMP = Z(J4 + 1) / Z(J4 - 2)
        Z(J4) = Z(J4 - 1) * TEMP
        D = D * TEMP
      else
        Z(J4) = Z(J4 + 1) * (Z(J4 - 1) / Z(J4 - 2))
        D = Z(J4 + 1) * (D / Z(J4 - 2))
      end if
      DMIN = MIN(DMIN, D)
      EMIN = MIN(EMIN, Z(J4))
    end do
  else
    do J4 = 4 * I0, 4 * (N0 - 3), 4
      Z(J4 - 3) = D + Z(J4)
      if (Z(J4 - 3) == ZERO) then
        Z(J4 - 1) = ZERO
        D = Z(J4 + 2)
        DMIN = D
        EMIN = ZERO
      else if (SAFMIN * Z(J4 + 2) < Z(J4 - 3) .and. SAFMIN * Z(J4 - 3) < Z(J4 + 2)) then
        TEMP = Z(J4 + 2) / Z(J4 - 3)
        Z(J4 - 1) = Z(J4) * TEMP
        D = D * TEMP
      else
        Z(J4 - 1) = Z(J4 + 2) * (Z(J4) / Z(J4 - 3))
        D = Z(J4 + 2) * (D / Z(J4 - 3))
      end if
      DMIN = MIN(DMIN, D)
      EMIN = MIN(EMIN, Z(J4 - 1))
    end do
  end if
!
! Unroll last two steps.
!
  DNM2 = D
  DMIN2 = DMIN
  J4 = 4 * (N0 - 2) - PP
  J4P2 = J4 + 2 * PP - 1
  Z(J4 - 2) = DNM2 + Z(J4P2)
  if (Z(J4 - 2) == ZERO) then
    Z(J4) = ZERO
    DNM1 = Z(J4P2 + 2)
    DMIN = DNM1
    EMIN = ZERO
  else if (SAFMIN * Z(J4P2 + 2) < Z(J4 - 2) .and. SAFMIN * Z(J4 - 2) < Z(J4P2 + 2)) then
    TEMP = Z(J4P2 + 2) / Z(J4 - 2)
    Z(J4) = Z(J4P2) * TEMP
    DNM1 = DNM2 * TEMP
  else
    Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
    DNM1 = Z(J4P2 + 2) * (DNM2 / Z(J4 - 2))
  end if
  DMIN = MIN(DMIN, DNM1)
!
  DMIN1 = DMIN
  J4 = J4 + 4
  J4P2 = J4 + 2 * PP - 1
  Z(J4 - 2) = DNM1 + Z(J4P2)
  if (Z(J4 - 2) == ZERO) then
    Z(J4) = ZERO
    DN = Z(J4P2 + 2)
    DMIN = DN
    EMIN = ZERO
  else if (SAFMIN * Z(J4P2 + 2) < Z(J4 - 2) .and. SAFMIN * Z(J4 - 2) < Z(J4P2 + 2)) then
    TEMP = Z(J4P2 + 2) / Z(J4 - 2)
    Z(J4) = Z(J4P2) * TEMP
    DN = DNM1 * TEMP
  else
    Z(J4) = Z(J4P2 + 2) * (Z(J4P2) / Z(J4 - 2))
    DN = Z(J4P2 + 2) * (DNM1 / Z(J4 - 2))
  end if
  DMIN = MIN(DMIN, DN)
!
  Z(J4 + 2) = DN
  Z(4 * N0 - PP) = EMIN
  return
!
! end of mobbrmsd_SLASQ6
!
end

