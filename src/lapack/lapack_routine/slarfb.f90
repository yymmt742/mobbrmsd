!> \brief \b SLARFB applies a block reflector or its transpose to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARFB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
!                          T, LDT, C, LDC, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               C( LDC, * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARFB applies a real block reflector H or its transpose H**T to a
!> real m by n matrix C, from either the left or the right.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply H or H**T from the Left
!>          = 'R': apply H or H**T from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply H (No transpose)
!>          = 'T': apply H**T (Transpose)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Indicates how H is formed from a product of elementary
!>          reflectors
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Indicates how the vectors which define the elementary
!>          reflectors are stored:
!>          = 'C': Columnwise
!>          = 'R': Rowwise
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the matrix T (= the number of elementary
!>          reflectors whose product defines the block reflector).
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is REAL array, dimension
!>                                (LDV,K) if STOREV = 'C'
!>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!>          The matrix V. See Further Details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!>          if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is REAL array, dimension (LDT,K)
!>          The triangular k by k matrix T in the representation of the
!>          block reflector.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is REAL array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LDWORK,K)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.
!>          If SIDE = 'L', LDWORK >= max(1,N);
!>          if SIDE = 'R', LDWORK >= max(1,M).
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
!> \date June 2013
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
!>  k = 3. The elements equal to 1 are not stored; the corresponding
!>  array elements are modified but restored on exit. The rest of the
!>  array is not used.
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
pure subroutine SLARFB(SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
     &                   T, LDT, C, LDC, WORK, LDWORK)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2013
!
!     .. Scalar Arguments ..
  character, intent(in) :: DIRECT, SIDE, STOREV, TRANS
  integer, intent(in)   :: K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout) :: C(LDC, *)
  real(RK), intent(in)    :: T(LDT, *), V(LDV, *)
  real(RK), intent(out)   :: WORK(LDWORK, *)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
  character :: TRANST
  integer   :: I, J
!
!     .. Parameters ..
! real(RK), parameter :: ONE = 1.0E+0
!..
! interface
! .. External Functions ..
!   include 'lsame.h'
! .. External Subroutines ..
!   include 'scopy.h'
!   include 'sgemm.h'
!   include 'strmm.h'
! end interface
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
  if (M <= 0 .or. N <= 0) return
!
  if (LSAME(TRANS, 'N')) then
    TRANST = 'T'
  else
    TRANST = 'N'
  end if
!
  if (LSAME(STOREV, 'C')) then
!
    if (LSAME(DIRECT, 'F')) then
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
      if (LSAME(SIDE, 'L')) then
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!              W := C1**T
!
        do J = 1, K
          call SCOPY(N, C(J, 1), LDC, WORK(1, J), 1)
        end do
!
!              W := W * V1
!
        call STRMM('Right', 'Lower', 'No transpose', 'Unit', N, &
&                     K, ONE, V, LDV, WORK, LDWORK)
        if (M > K) then
!
!                 W := W + C2**T * V2
!
          call SGEMM('Transpose', 'No transpose', N, K, M - K,&
&                        ONE, C(K + 1, 1), LDC, V(K + 1, 1), LDV, &
&                        ONE, WORK, LDWORK)
        end if
!
!              W := W * T**T  or  W * T
!
        call STRMM('Right', 'Upper', TRANST, 'Non-unit', N, K,&
&                     ONE, T, LDT, WORK, LDWORK)
!
!              C := C - V * W**T
!
        if (M > K) then
!
!                 C2 := C2 - V2 * W**T
!
          call SGEMM('No transpose', 'Transpose', M - K, N, K, &
&                        -ONE, V(K + 1, 1), LDV, WORK, LDWORK, ONE, &
&                        C(K + 1, 1), LDC)
        end if
!
!              W := W * V1**T
!
        call STRMM('Right', 'Lower', 'Transpose', 'Unit', N, K, &
&                     ONE, V, LDV, WORK, LDWORK)
!
!              C1 := C1 - W**T
!
        do J = 1, K
          do I = 1, N
            C(J, I) = C(J, I) - WORK(I, J)
          end do
        end do
!
      else if (LSAME(SIDE, 'R')) then
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
        do J = 1, K
          call SCOPY(M, C(1, J), 1, WORK(1, J), 1)
        end do
!
!              W := W * V1
!
        call STRMM('Right', 'Lower', 'No transpose', 'Unit', M, &
&                     K, ONE, V, LDV, WORK, LDWORK)
        if (N > K) then
!
!                 W := W + C2 * V2
!
          call SGEMM('No transpose', 'No transpose', M, K, N - K, &
&                        ONE, C(1, K + 1), LDC, V(K + 1, 1), LDV, &
&                        ONE, WORK, LDWORK)
        end if
!
!              W := W * T  or  W * T**T
!
        call STRMM('Right', 'Upper', TRANS, 'Non-unit', M, K, &
&                     ONE, T, LDT, WORK, LDWORK)
!
!              C := C - W * V**T
!
        if (N > K) then
!
!                 C2 := C2 - W * V2**T
!
          call SGEMM('No transpose', 'Transpose', M, N - K, K, &
&                        -ONE, WORK, LDWORK, V(K + 1, 1), LDV, ONE, &
&                        C(1, K + 1), LDC)
        end if
!
!              W := W * V1**T
!
        call STRMM('Right', 'Lower', 'Transpose', 'Unit', M, K, &
&                     ONE, V, LDV, WORK, LDWORK)
!
!              C1 := C1 - W
!
        do J = 1, K
          do I = 1, M
            C(I, J) = C(I, J) - WORK(I, J)
          end do
        end do
      end if
!
    else
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
      if (LSAME(SIDE, 'L')) then
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!              W := C2**T
!
        do J = 1, K
          call SCOPY(N, C(M - K + J, 1), LDC, WORK(1, J), 1)
        end do
!
!              W := W * V2
!
        call STRMM('Right', 'Upper', 'No transpose', 'Unit', N, &
&                     K, ONE, V(M - K + 1, 1), LDV, WORK, LDWORK)
        if (M > K) then
!
!                 W := W + C1**T * V1
!
          call SGEMM('Transpose', 'No transpose', N, K, M - K, &
&                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK)
        end if
!
!              W := W * T**T  or  W * T
!
        call STRMM('Right', 'Lower', TRANST, 'Non-unit', N, K, &
&                     ONE, T, LDT, WORK, LDWORK)
!
!              C := C - V * W**T
!
        if (M > K) then
!
!                 C1 := C1 - V1 * W**T
!
          call SGEMM('No transpose', 'Transpose', M - K, N, K, &
&                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC)
        end if
!
!              W := W * V2**T
!
        call STRMM('Right', 'Upper', 'Transpose', 'Unit', N, K, &
&                     ONE, V(M - K + 1, 1), LDV, WORK, LDWORK)
!
!              C2 := C2 - W**T
!
        do J = 1, K
          do I = 1, N
            C(M - K + J, I) = C(M - K + J, I) - WORK(I, J)
          end do
        end do
!
      else if (LSAME(SIDE, 'R')) then
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
        do J = 1, K
          call SCOPY(M, C(1, N - K + J), 1, WORK(1, J), 1)
        end do
!
!              W := W * V2
!
        call STRMM('Right', 'Upper', 'No transpose', 'Unit', M, &
&                     K, ONE, V(N - K + 1, 1), LDV, WORK, LDWORK)
        if (N > K) then
!
!                 W := W + C1 * V1
!
          call SGEMM('No transpose', 'No transpose', M, K, N - K, &
&                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK)
        end if
!
!              W := W * T  or  W * T**T
!
        call STRMM('Right', 'Lower', TRANS, 'Non-unit', M, K, &
&                     ONE, T, LDT, WORK, LDWORK)
!
!              C := C - W * V**T
!
        if (N > K) then
!
!                 C1 := C1 - W * V1**T
!
          call SGEMM('No transpose', 'Transpose', M, N - K, K, &
&                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC)
        end if
!
!              W := W * V2**T
!
        call STRMM('Right', 'Upper', 'Transpose', 'Unit', M, K, &
&                     ONE, V(N - K + 1, 1), LDV, WORK, LDWORK)
!
!              C2 := C2 - W
!
        do J = 1, K
          do I = 1, M
            C(I, N - K + J) = C(I, N - K + J) - WORK(I, J)
          end do
        end do
      end if
    end if
!
  else if (LSAME(STOREV, 'R')) then
!
    if (LSAME(DIRECT, 'F')) then
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
      if (LSAME(SIDE, 'L')) then
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!              W := C1**T
!
        do J = 1, K
          call SCOPY(N, C(J, 1), LDC, WORK(1, J), 1)
        end do
!
!              W := W * V1**T
!
        call STRMM('Right', 'Upper', 'Transpose', 'Unit', N, K, &
&                     ONE, V, LDV, WORK, LDWORK)
        if (M > K) then
!
!                 W := W + C2**T * V2**T
!
          call SGEMM('Transpose', 'Transpose', N, K, M - K, ONE, &
&                        C(K + 1, 1), LDC, V(1, K + 1), LDV, ONE, &
&                        WORK, LDWORK)
        end if
!
!              W := W * T**T  or  W * T
!
        call STRMM('Right', 'Upper', TRANST, 'Non-unit', N, K, &
&                     ONE, T, LDT, WORK, LDWORK)
!
!              C := C - V**T * W**T
!
        if (M > K) then
!
!                 C2 := C2 - V2**T * W**T
!
          call SGEMM('Transpose', 'Transpose', M - K, N, K, -ONE, &
&                        V(1, K + 1), LDV, WORK, LDWORK, ONE, &
&                        C(K + 1, 1), LDC)
        end if
!
!              W := W * V1
!
        call STRMM('Right', 'Upper', 'No transpose', 'Unit', N, &
&                     K, ONE, V, LDV, WORK, LDWORK)
!
!              C1 := C1 - W**T
!
        do J = 1, K
          do I = 1, N
            C(J, I) = C(J, I) - WORK(I, J)
          end do
        end do
!
      else if (LSAME(SIDE, 'R')) then
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C1
!
        do J = 1, K
          call SCOPY(M, C(1, J), 1, WORK(1, J), 1)
        end do
!
!              W := W * V1**T
!
        call STRMM('Right', 'Upper', 'Transpose', 'Unit', M, K, &
&                     ONE, V, LDV, WORK, LDWORK)
        if (N > K) then
!
!                 W := W + C2 * V2**T
!
          call SGEMM('No transpose', 'Transpose', M, K, N - K, &
&                        ONE, C(1, K + 1), LDC, V(1, K + 1), LDV, &
&                        ONE, WORK, LDWORK)
        end if
!
!              W := W * T  or  W * T**T
!
        call STRMM('Right', 'Upper', TRANS, 'Non-unit', M, K, &
&                     ONE, T, LDT, WORK, LDWORK)
!
!              C := C - W * V
!
        if (N > K) then
!
!                 C2 := C2 - W * V2
!
          call SGEMM('No transpose', 'No transpose', M, N - K, K, &
&                        -ONE, WORK, LDWORK, V(1, K + 1), LDV, ONE, &
&                        C(1, K + 1), LDC)
        end if
!
!              W := W * V1
!
        call STRMM('Right', 'Upper', 'No transpose', 'Unit', M, &
&                     K, ONE, V, LDV, WORK, LDWORK)
!
!              C1 := C1 - W
!
        do J = 1, K
          do I = 1, M
            C(I, J) = C(I, J) - WORK(I, J)
          end do
        end do
!
      end if
!
    else
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
      if (LSAME(SIDE, 'L')) then
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!              W := C2**T
!
        do J = 1, K
          call SCOPY(N, C(M - K + J, 1), LDC, WORK(1, J), 1)
        end do
!
!              W := W * V2**T
!
        call STRMM('Right', 'Lower', 'Transpose', 'Unit', N, K, &
&                     ONE, V(1, M - K + 1), LDV, WORK, LDWORK)
        if (M > K) then
!
!                 W := W + C1**T * V1**T
!
          call SGEMM('Transpose', 'Transpose', N, K, M - K, ONE, &
&                        C, LDC, V, LDV, ONE, WORK, LDWORK)
        end if
!
!              W := W * T**T  or  W * T
!
        call STRMM('Right', 'Lower', TRANST, 'Non-unit', N, K, &
&                     ONE, T, LDT, WORK, LDWORK)
!
!              C := C - V**T * W**T
!
        if (M > K) then
!
!                 C1 := C1 - V1**T * W**T
!
          call SGEMM('Transpose', 'Transpose', M - K, N, K, -ONE, &
              &       V, LDV, WORK, LDWORK, ONE, C, LDC)
        end if
!
!              W := W * V2
!
        call STRMM('Right', 'Lower', 'No transpose', 'Unit', N, &
            &      K, ONE, V(1, M - K + 1), LDV, WORK, LDWORK)
!
!              C2 := C2 - W**T
!
        do J = 1, K
          do I = 1, N
            C(M - K + J, I) = C(M - K + J, I) - WORK(I, J)
          end do
        end do
!
      else if (LSAME(SIDE, 'R')) then
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C2
!
        do J = 1, K
          call SCOPY(M, C(1, N - K + J), 1, WORK(1, J), 1)
        end do
!
!              W := W * V2**T
!
        call STRMM('Right', 'Lower', 'Transpose', 'Unit', M, K, &
            &      ONE, V(1, N - K + 1), LDV, WORK, LDWORK)
        if (N > K) then
!
!                 W := W + C1 * V1**T
!
          call SGEMM('No transpose', 'Transpose', M, K, N - K, &
              &      ONE, C, LDC, V, LDV, ONE, WORK, LDWORK)
        end if
!
!              W := W * T  or  W * T**T
!
        call STRMM('Right', 'Lower', TRANS, 'Non-unit', M, K, &
            &      ONE, T, LDT, WORK, LDWORK)
!
!              C := C - W * V
!
        if (N > K) then
!
!                 C1 := C1 - W * V1
!
          call SGEMM('No transpose', 'No transpose', M, N - K, K, &
              &      -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC)
        end if
!
!              W := W * V2
!
        call STRMM('Right', 'Lower', 'No transpose', 'Unit', M, &
            &      K, ONE, V(1, N - K + 1), LDV, WORK, LDWORK)
!
!              C1 := C1 - W
!
        do J = 1, K
          do I = 1, M
            C(I, N - K + J) = C(I, N - K + J) - WORK(I, J)
          end do
        end do
      end if
    end if
  end if
!
  return
!
!     End of SLARFB
!
end
