!| mobbrmsd_DLARFB applies a real block reflector H or its transpose H**T to a
!  real m by n matrix C, from either the left or the right.
!
!   The shape of the matrix V and the storage of the vectors which define
!   the H(i) is best illustrated by the following example with n = 5 and
!   k = 3. The elements equal to 1 are not stored; the corresponding
!   array elements are modified but restored on exit. The rest of the
!   array is not used.
!
!   DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!                V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                    ( v1  1    )                     (     1 v2 v2 v2 )
!                    ( v1 v2  1 )                     (        1 v3 v3 )
!                    ( v1 v2 v3 )
!                    ( v1 v2 v3 )
!
!   DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!                V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                    ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                    (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                    (     1 v3 )
!                    (        1 )
!
!  Reference DLARFB is provided by [netlib.org](http://www.netlib.org/lapack/).
!
!  -- LAPACK auxiliary routine --
!
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
pure subroutine mobbrmsd_DLARFB(SIDE, TRANS, DIRECT, STOREV, M, N, K, V, &
    &                           LDV, T, LDT, C, LDC, WORK, LDWORK)
  implicit none
  character, intent(in)    :: SIDE
!!          = 'L': apply H or H**T from the Left
!!
!!          = 'R': apply H or H**T from the Right
!!
  character, intent(in)    :: TRANS
!!          = 'N': apply H (No transpose)
!!
!!          = 'T': apply H**T (Transpose)
!!
  character, intent(in)    :: DIRECT
!!          Indicates how H is formed from a product of elementary
!!          reflectors
!!
!!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!!
!!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!!
  character, intent(in)    :: STOREV
!!          Indicates how the vectors which define the elementary
!!          reflectors are stored:
!!
!!          = 'C': Columnwise
!!
!!          = 'R': Rowwise
!!
  integer, intent(in)      :: M
!!          The number of rows of the matrix C.
!!
  integer, intent(in)      :: N
!!          The number of columns of the matrix C.
!!
  integer, intent(in)      :: K
!!          The order of the matrix T (= the number of elementary
!!          reflectors whose product defines the block reflector).
!!
!!          If SIDE = 'L', M >= K >= 0.
!!
!!          If SIDE = 'R', N >= K >= 0.
!!
  integer, intent(in)      :: LDV
!!          The leading dimension of the array V.
!!
!!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M).
!!
!!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N).
!!
!!          if STOREV = 'R', LDV >= K.
!!
  real(RK), intent(in)     :: V(LDV, *)
!!          V is real(RK)           :: array, dimension
!!
!!          The matrix V. See Further Details.
!!
!!                                (LDV,K) if STOREV = 'C'
!!
!!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!!
!!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!!
  integer, intent(in)      :: LDT
!!          The leading dimension of the array T. LDT >= K.
!!
  real(RK), intent(in)     :: T(LDT, *)
!!          T is real(RK)           :: array, dimension (LDT,K)
!!
!!          The triangular k by k matrix T in the representation of the
!!          block reflector.
!!
  integer, intent(in)      :: LDC
!!          The leading dimension of the array C. LDC >= max(1,M).
!!
  real(RK), intent(inout)  :: C(LDC, *)
!!          C is real(RK)           :: array, dimension (LDC,N)
!!
!!          On entry, the m by n matrix C.
!!
!!          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.
!!
  integer, intent(in)      :: LDWORK
!!          The leading dimension of the array WORK.
!!
!!          If SIDE = 'L', LDWORK >= max(1,N);
!!
!!          if SIDE = 'R', LDWORK >= max(1,M).
!!
  real(RK), intent(out)    :: WORK(LDWORK, *)
!!          WORK is real(RK)           :: array, dimension (LDWORK,K)
!!
  character               :: TRANST
  integer                 :: I, J
! real(RK), parameter     :: ONE = 1.0_RK
! interface
!     .. External Subroutines ..
!   include 'dcopy.h'
!   include 'dgemm.h'
!   include 'dtrmm.h'
!     .. External Functions ..
!   include 'lsame.h'
! end interface
!
! Quick return if possible
!
  if (M <= 0 .or. N <= 0) return
!
  if (mobbrmsd_LSAME(TRANS, 'N')) then
    TRANST = 'T'
  else
    TRANST = 'N'
  end if
!
  if (mobbrmsd_LSAME(STOREV, 'C')) then
!
    if (mobbrmsd_LSAME(DIRECT, 'F')) then
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
      if (mobbrmsd_LSAME(SIDE, 'L')) then
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!              W := C1**T
!
        do J = 1, K
          call mobbrmsd_DCOPY(N, C(J, 1), LDC, WORK(1, J), 1)
        end do
!
!            W := W * V1
!
        call mobbrmsd_DTRMM('Right', 'Lower', 'No transpose', 'Unit', N,&
       &            K, ONE, V, LDV, WORK, LDWORK)
        if (M > K) then
!
!               W := W + C2**T * V2
!
          call mobbrmsd_DGEMM('Transpose', 'No transpose', N, K, M - K,&
         &            ONE, C(K + 1, 1), LDC, V(K + 1, 1), LDV,&
         &            ONE, WORK, LDWORK)
        end if
!
!            W := W * T**T  or  W * T
!
        call mobbrmsd_DTRMM('Right', 'Upper', TRANST, 'Non-unit', N, K,&
       &            ONE, T, LDT, WORK, LDWORK)
!
!            C := C - V * W**T
!
        if (M > K) then
!
!               C2 := C2 - V2 * W**T
!
          call mobbrmsd_DGEMM('No transpose', 'Transpose', M - K, N, K,&
         &            -ONE, V(K + 1, 1), LDV, WORK, LDWORK, ONE,&
         &            C(K + 1, 1), LDC)
        end if
!
!            W := W * V1**T
!
        call mobbrmsd_DTRMM('Right', 'Lower', 'Transpose', 'Unit', N, K,&
       &            ONE, V, LDV, WORK, LDWORK)
!
!            C1 := C1 - W**T
!
!             do 30 J = 1, K
!               do 20 I = 1, N
!                 C(J, I) = C(J, I) - WORK(I, J)
!20             continue
!30           continue
        do concurrent(J=1:K, I=1:N)
          C(J, I) = C(J, I) - WORK(I, J)
        end do
!
      else if (mobbrmsd_LSAME(SIDE, 'R')) then
!
!        Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!        W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!        W := C1
!
        do J = 1, K
          call mobbrmsd_DCOPY(M, C(1, J), 1, WORK(1, J), 1)
        end do
!
!        W := W * V1
!
        call mobbrmsd_DTRMM('Right', 'Lower', 'No transpose', 'Unit', M,&
       &            K, ONE, V, LDV, WORK, LDWORK)
        if (N > K) then
!
!             W := W + C2 * V2
!
          call mobbrmsd_DGEMM('No transpose', 'No transpose', M, K, N - K,&
         &            ONE, C(1, K + 1), LDC, V(K + 1, 1), LDV,&
         &            ONE, WORK, LDWORK)
        end if
!
!          W := W * T  or  W * T**T
!
        call mobbrmsd_DTRMM('Right', 'Upper', TRANS, 'Non-unit', M, K,&
       &            ONE, T, LDT, WORK, LDWORK)
!
!          C := C - W * V**T
!
        if (N > K) then
!
!             C2 := C2 - W * V2**T
!
          call mobbrmsd_DGEMM('No transpose', 'Transpose', M, N - K, K,&
         &            -ONE, WORK, LDWORK, V(K + 1, 1), LDV, ONE,&
         &            C(1, K + 1), LDC)
        end if
!
!          W := W * V1**T
!
        call mobbrmsd_DTRMM('Right', 'Lower', 'Transpose', 'Unit', M, K,&
       &            ONE, V, LDV, WORK, LDWORK)
!
!          C1 := C1 - W
!
        do concurrent(J=1:K, I=1:M)
          C(I, J) = C(I, J) - WORK(I, J)
        end do
!           do 60 J = 1, K
!             do 50 I = 1, M
!               C(I, J) = C(I, J) - WORK(I, J)
!50           continue
!60         continue
      end if
!
    else
!
!       Let  V =  ( V1 )
!                 ( V2 )    (last K rows)
!       where  V2  is unit upper triangular.
!
      if (mobbrmsd_LSAME(SIDE, 'L')) then
!
!          Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                ( C2 )
!
!          W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!          W := C2**T
!
        do J = 1, K
          call mobbrmsd_DCOPY(N, C(M - K + J, 1), LDC, WORK(1, J), 1)
        end do
!
!          W := W * V2
!
        call mobbrmsd_DTRMM('Right', 'Upper', 'No transpose', 'Unit', N,&
       &            K, ONE, V(M - K + 1, 1), LDV, WORK, LDWORK)
        if (M > K) then
!
!             W := W + C1**T * V1
!
          call mobbrmsd_DGEMM('Transpose', 'No transpose', N, K, M - K,&
         &            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK)
        end if
!
!          W := W * T**T  or  W * T
!
        call mobbrmsd_DTRMM('Right', 'Lower', TRANST, 'Non-unit', N, K,&
       &            ONE, T, LDT, WORK, LDWORK)
!
!          C := C - V * W**T
!
        if (M > K) then
!
!             C1 := C1 - V1 * W**T
!
          call mobbrmsd_DGEMM('No transpose', 'Transpose', M - K, N, K,&
         &            -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC)
        end if
!
!          W := W * V2**T
!
        call mobbrmsd_DTRMM('Right', 'Upper', 'Transpose', 'Unit', N, K,&
       &            ONE, V(M - K + 1, 1), LDV, WORK, LDWORK)
!
!          C2 := C2 - W**T
!
        do concurrent(J=1:K, I=1:N)
          C(M - K + J, I) = C(M - K + J, I) - WORK(I, J)
        end do
!           do 90 J = 1, K
!             do 80 I = 1, N
!               C(M - K + J, I) = C(M - K + J, I) - WORK(I, J)
!80           continue
!90         continue
!
      else if (mobbrmsd_LSAME(SIDE, 'R')) then
!
!           Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!           W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!           W := C2
!
        do J = 1, K
          call mobbrmsd_DCOPY(M, C(1, N - K + J), 1, WORK(1, J), 1)
        end do
!
!          W := W * V2
!
        call mobbrmsd_DTRMM('Right', 'Upper', 'No transpose', 'Unit', M,&
       &            K, ONE, V(N - K + 1, 1), LDV, WORK, LDWORK)
        if (N > K) then
!
! W := W + C1 * V1
!
          call mobbrmsd_DGEMM('No transpose', 'No transpose', M, K, N - K,&
         &            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK)
        end if
!
!         W := W * T  or  W * T**T
!
        call mobbrmsd_DTRMM('Right', 'Lower', TRANS, 'Non-unit', M, K,&
       &            ONE, T, LDT, WORK, LDWORK)
!
!           C := C - W * V**T
!
        if (N > K) then
!
!             C1 := C1 - W * V1**T
!
          call mobbrmsd_DGEMM('No transpose', 'Transpose', M, N - K, K,&
         &            -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC)
        end if
!
!         W := W * V2**T
!
        call mobbrmsd_DTRMM('Right', 'Upper', 'Transpose', 'Unit', M, K,&
       &            ONE, V(N - K + 1, 1), LDV, WORK, LDWORK)
!
!        C2 := C2 - W
!
        do concurrent(J=1:K, I=1:M)
          C(I, N - K + J) = C(I, N - K + J) - WORK(I, J)
        end do
!             do 120 J = 1, K
!               do 110 I = 1, M
!                 C(I, N - K + J) = C(I, N - K + J) - WORK(I, J)
!110            continue
!120          continue
      end if
    end if
!
  else if (mobbrmsd_LSAME(STOREV, 'R')) then
!
    if (mobbrmsd_LSAME(DIRECT, 'F')) then
!
!       Let  V =  ( V1  V2 )    (V1: first K columns)
!       where  V1  is unit upper triangular.
!
      if (mobbrmsd_LSAME(SIDE, 'L')) then
!
!          Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                ( C2 )
!
!          W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!          W := C1**T
!
        do J = 1, K
          call mobbrmsd_DCOPY(N, C(J, 1), LDC, WORK(1, J), 1)
        end do
!
!          W := W * V1**T
!
        call mobbrmsd_DTRMM('Right', 'Upper', 'Transpose', 'Unit', N, K,&
       &            ONE, V, LDV, WORK, LDWORK)
        if (M > K) then
!
!             W := W + C2**T * V2**T
!
          call mobbrmsd_DGEMM('Transpose', 'Transpose', N, K, M - K, ONE,&
         &            C(K + 1, 1), LDC, V(1, K + 1), LDV, ONE,&
         &            WORK, LDWORK)
        end if
!
!          W := W * T**T  or  W * T
!
        call mobbrmsd_DTRMM('Right', 'Upper', TRANST, 'Non-unit', N, K,&
       &            ONE, T, LDT, WORK, LDWORK)
!
!          C := C - V**T * W**T
!
        if (M > K) then
!
!             C2 := C2 - V2**T * W**T
!
          call mobbrmsd_DGEMM('Transpose', 'Transpose', M - K, N, K, -ONE,&
         &            V(1, K + 1), LDV, WORK, LDWORK, ONE,&
         &            C(K + 1, 1), LDC)
        end if
!
!          W := W * V1
!
        call mobbrmsd_DTRMM('Right', 'Upper', 'No transpose', 'Unit', N,&
       &            K, ONE, V, LDV, WORK, LDWORK)
!
!          C1 := C1 - W**T
!
        do concurrent(J=1:K, I=1:N)
          C(J, I) = C(J, I) - WORK(I, J)
        end do
!             do 150 J = 1, K
!               do 140 I = 1, N
!                 C(J, I) = C(J, I) - WORK(I, J)
!140            continue
!150          continue
!
      else if (mobbrmsd_LSAME(SIDE, 'R')) then
!
!          Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!          W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!          W := C1
!
        do J = 1, K
          call mobbrmsd_DCOPY(M, C(1, J), 1, WORK(1, J), 1)
        end do
!
!          W := W * V1**T
!
        call mobbrmsd_DTRMM('Right', 'Upper', 'Transpose', 'Unit', M, K,&
       &            ONE, V, LDV, WORK, LDWORK)
        if (N > K) then
!
!             W := W + C2 * V2**T
!
          call mobbrmsd_DGEMM('No transpose', 'Transpose', M, K, N - K,&
         &            ONE, C(1, K + 1), LDC, V(1, K + 1), LDV,&
         &            ONE, WORK, LDWORK)
        end if
!
!          W := W * T  or  W * T**T
!
        call mobbrmsd_DTRMM('Right', 'Upper', TRANS, 'Non-unit', M, K,&
       &            ONE, T, LDT, WORK, LDWORK)
!
!          C := C - W * V
!
        if (N > K) then
!
!             C2 := C2 - W * V2
!
          call mobbrmsd_DGEMM('No transpose', 'No transpose', M, N - K, K,&
         &            -ONE, WORK, LDWORK, V(1, K + 1), LDV, ONE,&
         &            C(1, K + 1), LDC)
        end if
!
!          W := W * V1
!
        call mobbrmsd_DTRMM('Right', 'Upper', 'No transpose', 'Unit', M,&
       &            K, ONE, V, LDV, WORK, LDWORK)
!
!          C1 := C1 - W
!
        do concurrent(I=1:M, J=1:K)
          C(I, J) = C(I, J) - WORK(I, J)
        end do
!         do 180 J = 1, K
!           do 170 I = 1, M
!             C(I, J) = C(I, J) - WORK(I, J)
!170        continue
!180      continue
!
      end if
!
    else
!
!       Let  V =  ( V1  V2 )    (V2: last K columns)
!       where  V2  is unit lower triangular.
!
      if (mobbrmsd_LSAME(SIDE, 'L')) then
!
!          Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                ( C2 )
!
!          W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!          W := C2**T
!
        do J = 1, K
          call mobbrmsd_DCOPY(N, C(M - K + J, 1), LDC, WORK(1, J), 1)
        end do
!
!          W := W * V2**T
!
        call mobbrmsd_DTRMM('Right', 'Lower', 'Transpose', 'Unit', N, K,&
       &            ONE, V(1, M - K + 1), LDV, WORK, LDWORK)
        if (M > K) then
!
!             W := W + C1**T * V1**T
!
          call mobbrmsd_DGEMM('Transpose', 'Transpose', N, K, M - K, ONE,&
         &            C, LDC, V, LDV, ONE, WORK, LDWORK)
        end if
!
!          W := W * T**T  or  W * T
!
        call mobbrmsd_DTRMM('Right', 'Lower', TRANST, 'Non-unit', N, K,&
       &            ONE, T, LDT, WORK, LDWORK)
!
!          C := C - V**T * W**T
!
        if (M > K) then
!
!             C1 := C1 - V1**T * W**T
!
          call mobbrmsd_DGEMM('Transpose', 'Transpose', M - K, N, K, -ONE,&
         &            V, LDV, WORK, LDWORK, ONE, C, LDC)
        end if
!
!          W := W * V2
!
        call mobbrmsd_DTRMM('Right', 'Lower', 'No transpose', 'Unit', N,&
       &            K, ONE, V(1, M - K + 1), LDV, WORK, LDWORK)
!
!          C2 := C2 - W**T
!
        do concurrent(J=1:K, I=1:N)
          C(M - K + J, I) = C(M - K + J, I) - WORK(I, J)
        end do
!             do 210 J = 1, K
!               do 200 I = 1, N
!                 C(M - K + J, I) = C(M - K + J, I) - WORK(I, J)
!200            continue
!210          continue
!
      else if (mobbrmsd_LSAME(SIDE, 'R')) then
!
!          Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!          W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!          W := C2
!
        do J = 1, K
          call mobbrmsd_DCOPY(M, C(1, N - K + J), 1, WORK(1, J), 1)
        end do
!
!          W := W * V2**T
!
        call mobbrmsd_DTRMM('Right', 'Lower', 'Transpose', 'Unit', M, K,&
       &            ONE, V(1, N - K + 1), LDV, WORK, LDWORK)
        if (N > K) then
!
!             W := W + C1 * V1**T
!
          call mobbrmsd_DGEMM('No transpose', 'Transpose', M, K, N - K,&
         &            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK)
        end if
!
!          W := W * T  or  W * T**T
!
        call mobbrmsd_DTRMM('Right', 'Lower', TRANS, 'Non-unit', M, K,&
       &            ONE, T, LDT, WORK, LDWORK)
!
!          C := C - W * V
!
        if (N > K) then
!
!             C1 := C1 - W * V1
!
          call mobbrmsd_DGEMM('No transpose', 'No transpose', M, N - K, K,&
         &            -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC)
        end if
!
!          W := W * V2
!
        call mobbrmsd_DTRMM('Right', 'Lower', 'No transpose', 'Unit', M,&
       &            K, ONE, V(1, N - K + 1), LDV, WORK, LDWORK)
!
!          C1 := C1 - W
!
        do concurrent(J=1:K, I=1:M)
          C(I, N - K + J) = C(I, N - K + J) - WORK(I, J)
        end do
!             do 240 J = 1, K
!               do 230 I = 1, M
!                 C(I, N - K + J) = C(I, N - K + J) - WORK(I, J)
!230            continue
!240          continue
!
      end if
!
    end if
  end if
!
  return
!
!       End of mobbrmsd_DLARFB
!
end subroutine mobbrmsd_DLARFB

