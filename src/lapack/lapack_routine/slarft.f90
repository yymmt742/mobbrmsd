!> \brief \b SLARFT forms the triangular factor T of a block reflector H = I - vtvH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARFT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarft.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarft.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarft.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, STOREV
!       INTEGER            K, LDT, LDV, N
!       ..
!       .. Array Arguments ..
!       REAL               T( LDT, * ), TAU( * ), V( LDV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARFT forms the triangular factor T of a real block reflector H
!> of order n, which is defined as a product of k elementary reflectors.
!>
!> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!>
!> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!>
!> If STOREV = 'C', the vector which defines the elementary reflector
!> H(i) is stored in the i-th column of the array V, and
!>
!>    H  =  I - V * T * V**T
!>
!> If STOREV = 'R', the vector which defines the elementary reflector
!> H(i) is stored in the i-th row of the array V, and
!>
!>    H  =  I - V**T * T * V
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Specifies the order in which the elementary reflectors are
!>          multiplied to form the block reflector:
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Specifies how the vectors which define the elementary
!>          reflectors are stored (see also Further Details):
!>          = 'C': columnwise
!>          = 'R': rowwise
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the block reflector H. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the triangular factor T (= the number of
!>          elementary reflectors). K >= 1.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is REAL array, dimension
!>                               (LDV,K) if STOREV = 'C'
!>                               (LDV,N) if STOREV = 'R'
!>          The matrix V. See further details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is REAL array, dimension (LDT,K)
!>          The k by k triangular factor T of the block reflector.
!>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!>          lower triangular. The rest of the array is not used.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup realOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The elements equal to 1 are not stored.
!>
!>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!>
!>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!>                   ( v1  1    )                     (     1 v2 v2 v2 )
!>                   ( v1 v2  1 )                     (        1 v3 v3 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!>
!>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!>                   (     1 v3 )
!>                   (        1 )
!> \endverbatim
!>
!  =====================================================================
pure subroutine SLARFT(DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT)
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
  character, intent(in) :: DIRECT, STOREV
  integer, intent(in)   :: K, LDT, LDV, N
!..
!..Array Arguments..
  real(RK), intent(in)  :: TAU(*), V(LDV, *)
  real(RK), intent(out) :: T(LDT, *)
!..
!
!  =====================================================================
!
!..Local Scalars..
  integer :: I, J, PREVLASTV, LASTV
!..
!..Parameters..
! real(RK), parameter :: ZERO = 0.0E0
! real(RK), parameter :: ONE = 1.0E0
!..
! interface
! .. External Functions ..
!   include 'lsame.h'
!..external Subroutines..
!   include 'sgemv.h'
!   include 'strmv.h'
! end interface
!..
!..Executable Statements..
!
! Quick return if possible
!
  if (N == 0) return
!
  if (LSAME(DIRECT, 'F')) then
    PREVLASTV = N
    do I = 1, K
      PREVLASTV = MAX(I, PREVLASTV)
      if (TAU(I) == ZERO) then
!
! H(i) = I
!
        do J = 1, I
          T(J, I) = ZERO
        end do
      else
!
! general case
!
        if (LSAME(STOREV, 'C')) then
! Skip any trailing zeros.
          do LASTV = N, I + 1, -1
            if (V(LASTV, I) /= ZERO) exit
          end do
          do J = 1, I - 1
            T(J, I) = -TAU(I) * V(I, J)
          end do
          J = MIN(LASTV, PREVLASTV)
!
! T(1:i - 1, i): = -tau(i) * V(i:j, 1:i - 1)**T * V(i:j, i)
!
          call SGEMV('Transpose', J - I, I - 1, -TAU(I), &
          & V(I + 1, 1), LDV, V(I + 1, I), 1, ONE, &
          & T(1, I), 1)
        else
! Skip any trailing zeros.
          do LASTV = N, I + 1, -1
            if (V(I, LASTV) /= ZERO) exit
          end do
          do J = 1, I - 1
            T(J, I) = -TAU(I) * V(J, I)
          end do
          J = MIN(LASTV, PREVLASTV)
!
! T(1:i - 1, i): = -tau(i) * V(1:i - 1, i:j) * V(i, i:j)**T
!
          call SGEMV('No transpose', I - 1, J - I, -TAU(I), &
              &      V(1, I + 1), LDV, V(I, I + 1), LDV, ONE, T(1, I), 1)
        end if
!
! T(1:i - 1, i): = T(1:i - 1, 1:i - 1) * T(1:i - 1, i)
!
        call STRMV('Upper', 'No transpose', 'Non-unit', I - 1, T, &
            &      LDT, T(1, I), 1)
        T(I, I) = TAU(I)
        if (I > 1) then
          PREVLASTV = MAX(PREVLASTV, LASTV)
        else
          PREVLASTV = LASTV
        end if
      end if
    end do
  else
    PREVLASTV = 1
    do I = K, 1, -1
      if (TAU(I) == ZERO) then
!
! H(i) = I
!
        do J = I, K
          T(J, I) = ZERO
        end do
      else
!
! general case
!
        if (I < K) then
          if (LSAME(STOREV, 'C')) then
! Skip any leading zeros.
            do LASTV = 1, I - 1
              if (V(LASTV, I) /= ZERO) exit
            end do
            do J = I + 1, K
              T(J, I) = -TAU(I) * V(N - K + I, J)
            end do
            J = MAX(LASTV, PREVLASTV)
!
! T(i + 1:k, i) = -tau(i) * V(j:n - k + i, i + 1:k)**T * V(j:n - k + i, i)
!
            call SGEMV('Transpose', N - K + I - J, K - I, -TAU(I), &
                &      V(J, I + 1), LDV, V(J, I), 1, ONE, &
                &      T(I + 1, I), 1)
          else
! Skip any leading zeros.
            do LASTV = 1, I - 1
              if (V(I, LASTV) /= ZERO) exit
            end do
            do J = I + 1, K
              T(J, I) = -TAU(I) * V(J, N - K + I)
            end do
            J = MAX(LASTV, PREVLASTV)
!
! T(i + 1:k, i) = -tau(i) * V(i + 1:k, j:n - k + i) * V(i, j:n - k + i)**T
!
            call SGEMV('No transpose', K - I, N - K + I - J, &
                &      -TAU(I), V(I + 1, J), LDV, V(I, J), LDV, &
                &      ONE, T(I + 1, I), 1)
          end if
!
! T(i + 1:k, i): = T(i + 1:k, i + 1:k) * T(i + 1:k, i)
!
          call STRMV('Lower', 'No transpose', 'Non-unit', K - I, &
              &      T(I + 1, I + 1), LDT, T(I + 1, I), 1)
          if (I > 1) then
            PREVLASTV = MIN(PREVLASTV, LASTV)
          else
            PREVLASTV = LASTV
          end if
        end if
        T(I, I) = TAU(I)
      end if
    end do
  end if
  return
!
! end of SLARFT
!
end
