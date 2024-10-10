!| mobbrmsd_SCOMBSSQ adds two scaled sum of squares quantities, V1 := V1 + V2.
!  That is,
!
!     V1_scale**2 * V1_sumsq := V1_scale**2 * V1_sumsq
!                             + V2_scale**2 * V2_sumsq
!
!  Reference SCOMBSSQ is provided by [netlib](http://www.netlib.org/lapack/explore-html/).
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2018
!
pure subroutine mobbrmsd_SCOMBSSQ(V1, V2)
  implicit none
  real, intent(inout) :: V1(2)
!!  REAL array, dimension (2).
!!  The first scaled sum.
!!  V1(1) = V1_scale, V1(2) = V1_sumsq.
!!
  real, intent(in)    :: V2(2)
!!  REAL array, dimension (2).
!!  The second scaled sum.
!!  V2(1) = V2_scale, V2(2) = V2_sumsq.
!
  if (V1(1) >= V2(1)) then
    if (V1(1) /= ZERO) then
      V1(2) = V1(2) + (V2(1) / V1(1))**2 * V2(2)
    else
      V1(2) = V1(2) + V2(2)
    end if
  else
    V1(2) = V2(2) + (V1(1) / V2(1))**2 * V1(2)
    V1(1) = V2(1)
  end if
  return
!
! End of mobbrmsd_SCOMBSSQ
!
end

