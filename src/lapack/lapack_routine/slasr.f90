!> \brief \b SLASR applies a sequence of plane rotations to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, PIVOT, SIDE
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), C( * ), S( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASR applies a sequence of plane rotations to a real matrix A,
!> from either the left or the right.
!>
!> When SIDE = 'L', the transformation takes the form
!>
!>    A := P*A
!>
!> and when SIDE = 'R', the transformation takes the form
!>
!>    A := A*P**T
!>
!> where P is an orthogonal matrix consisting of a sequence of z plane
!> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
!> and P**T is the transpose of P.
!>
!> When DIRECT = 'F' (Forward sequence), then
!>
!>    P = P(z-1) * ... * P(2) * P(1)
!>
!> and when DIRECT = 'B' (Backward sequence), then
!>
!>    P = P(1) * P(2) * ... * P(z-1)
!>
!> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
!>
!>    R(k) = (  c(k)  s(k) )
!>         = ( -s(k)  c(k) ).
!>
!> When PIVOT = 'V' (Variable pivot), the rotation is performed
!> for the plane (k,k+1), i.e., P(k) has the form
!>
!>    P(k) = (  1                                            )
!>           (       ...                                     )
!>           (              1                                )
!>           (                   c(k)  s(k)                  )
!>           (                  -s(k)  c(k)                  )
!>           (                                1              )
!>           (                                     ...       )
!>           (                                            1  )
!>
!> where R(k) appears as a rank-2 modification to the identity matrix in
!> rows and columns k and k+1.
!>
!> When PIVOT = 'T' (Top pivot), the rotation is performed for the
!> plane (1,k+1), so P(k) has the form
!>
!>    P(k) = (  c(k)                    s(k)                 )
!>           (         1                                     )
!>           (              ...                              )
!>           (                     1                         )
!>           ( -s(k)                    c(k)                 )
!>           (                                 1             )
!>           (                                      ...      )
!>           (                                             1 )
!>
!> where R(k) appears in rows and columns 1 and k+1.
!>
!> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
!> performed for the plane (k,z), giving P(k) the form
!>
!>    P(k) = ( 1                                             )
!>           (      ...                                      )
!>           (             1                                 )
!>           (                  c(k)                    s(k) )
!>           (                         1                     )
!>           (                              ...              )
!>           (                                     1         )
!>           (                 -s(k)                    c(k) )
!>
!> where R(k) appears in rows and columns k and z.  The rotations are
!> performed without ever forming P(k) explicitly.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          Specifies whether the plane rotation matrix P is applied to
!>          A on the left or the right.
!>          = 'L':  Left, compute A := P*A
!>          = 'R':  Right, compute A:= A*P**T
!> \endverbatim
!>
!> \param[in] PIVOT
!> \verbatim
!>          PIVOT is CHARACTER*1
!>          Specifies the plane for which P(k) is a plane rotation
!>          matrix.
!>          = 'V':  Variable pivot, the plane (k,k+1)
!>          = 'T':  Top pivot, the plane (1,k+1)
!>          = 'B':  Bottom pivot, the plane (k,z)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Specifies whether P is a forward or backward sequence of
!>          plane rotations.
!>          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
!>          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  If m <= 1, an immediate
!>          return is effected.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  If n <= 1, an
!>          immediate return is effected.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension
!>                  (M-1) if SIDE = 'L'
!>                  (N-1) if SIDE = 'R'
!>          The cosines c(k) of the plane rotations.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension
!>                  (M-1) if SIDE = 'L'
!>                  (N-1) if SIDE = 'R'
!>          The sines s(k) of the plane rotations.  The 2-by-2 plane
!>          rotation part of the matrix P(k), R(k), has the form
!>          R(k) = (  c(k)  s(k) )
!>                 ( -s(k)  c(k) ).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The M-by-N matrix A.  On exit, A is overwritten by P*A if
!>          SIDE = 'R' or by A*P**T if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
pure subroutine SLASR(SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA)
  implicit none
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
  character, intent(in) :: DIRECT, PIVOT, SIDE
  integer, intent(in)   :: LDA, M, N
!..
!..Array Arguments..
  real(RK), intent(in)    :: C(*), S(*)
  real(RK), intent(inout) :: A(LDA, *)
!..
!
!  =====================================================================
!
!..Parameters..
! real, parameter :: ZERO = 0.0E0
! real, parameter :: ONE = 1.0E+0
!..
! interface
!..external Functions..
!   include 'lsame.h'
! end interface
!..Local Scalars..
  integer :: I, INFO, J
  real(RK) :: CTEMP, STEMP, TEMP
!..
!..intrinsic Functions..
  intrinsic :: MAX
!..
!..Executable Statements..
!
!Test the input parameters
!
  INFO = 0
  if (.not. (LSAME(SIDE, 'L') .or. LSAME(SIDE, 'R'))) then
    INFO = 1
  else if (.not. (LSAME(PIVOT, 'V') .or. LSAME(PIVOT, 'T') .or. LSAME(PIVOT, 'B'))) then
    INFO = 2
  else if (.not. (LSAME(DIRECT, 'F') .or. LSAME(DIRECT, 'B'))) then
    INFO = 3
  else if (M < 0) then
    INFO = 4
  else if (N < 0) then
    INFO = 5
  else if (LDA < MAX(1, M)) then
    INFO = 9
  end if
  if (INFO /= 0) then
!   call XERBLA('SLASR ', INFO)
    return
  end if
!
!Quick return if possible
!
  if ((M == 0) .or. (N == 0)) return
  if (LSAME(SIDE, 'L')) then
    !
    !Form P * A
    !
    if (LSAME(PIVOT, 'V')) then
      if (LSAME(DIRECT, 'F')) then
        do J = 1, M - 1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J + 1, I)
              A(J + 1, I) = CTEMP * TEMP - STEMP * A(J, I)
              A(J, I) = STEMP * TEMP + CTEMP * A(J, I)
            end do
          end if
        end do
      else if (LSAME(DIRECT, 'B')) then
        do J = M - 1, 1, -1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J + 1, I)
              A(J + 1, I) = CTEMP * TEMP - STEMP * A(J, I)
              A(J, I) = STEMP * TEMP + CTEMP * A(J, I)
            end do
          end if
        end do
      end if
    else if (LSAME(PIVOT, 'T')) then
      if (LSAME(DIRECT, 'F')) then
        do J = 2, M
          CTEMP = C(J - 1)
          STEMP = S(J - 1)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J, I)
              A(J, I) = CTEMP * TEMP - STEMP * A(1, I)
              A(1, I) = STEMP * TEMP + CTEMP * A(1, I)
            end do
          end if
        end do
      else if (LSAME(DIRECT, 'B')) then
        do J = M, 2, -1
          CTEMP = C(J - 1)
          STEMP = S(J - 1)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J, I)
              A(J, I) = CTEMP * TEMP - STEMP * A(1, I)
              A(1, I) = STEMP * TEMP + CTEMP * A(1, I)
            end do
          end if
        end do
      end if
    else if (LSAME(PIVOT, 'B')) then
      if (LSAME(DIRECT, 'F')) then
        do J = 1, M - 1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J, I)
              A(J, I) = STEMP * A(M, I) + CTEMP * TEMP
              A(M, I) = CTEMP * A(M, I) - STEMP * TEMP
            end do
          end if
        end do
      else if (LSAME(DIRECT, 'B')) then
        do J = M - 1, 1, -1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, N
              TEMP = A(J, I)
              A(J, I) = STEMP * A(M, I) + CTEMP * TEMP
              A(M, I) = CTEMP * A(M, I) - STEMP * TEMP
            end do
          end if
        end do
      end if
    end if
  else if (LSAME(SIDE, 'R')) then
    !
    !Form A * P**T
    !
    if (LSAME(PIVOT, 'V')) then
      if (LSAME(DIRECT, 'F')) then
        do J = 1, N - 1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J + 1)
              A(I, J + 1) = CTEMP * TEMP - STEMP * A(I, J)
              A(I, J) = STEMP * TEMP + CTEMP * A(I, J)
            end do
          end if
        end do
      else if (LSAME(DIRECT, 'B')) then
        do J = N - 1, 1, -1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J + 1)
              A(I, J + 1) = CTEMP * TEMP - STEMP * A(I, J)
              A(I, J) = STEMP * TEMP + CTEMP * A(I, J)
            end do
          end if
        end do
      end if
    else if (LSAME(PIVOT, 'T')) then
      if (LSAME(DIRECT, 'F')) then
        do J = 2, N
          CTEMP = C(J - 1)
          STEMP = S(J - 1)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J)
              A(I, J) = CTEMP * TEMP - STEMP * A(I, 1)
              A(I, 1) = STEMP * TEMP + CTEMP * A(I, 1)
            end do
          end if
        end do
      else if (LSAME(DIRECT, 'B')) then
        do J = N, 2, -1
          CTEMP = C(J - 1)
          STEMP = S(J - 1)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J)
              A(I, J) = CTEMP * TEMP - STEMP * A(I, 1)
              A(I, 1) = STEMP * TEMP + CTEMP * A(I, 1)
            end do
          end if
        end do
      end if
    else if (LSAME(PIVOT, 'B')) then
      if (LSAME(DIRECT, 'F')) then
        do J = 1, N - 1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J)
              A(I, J) = STEMP * A(I, N) + CTEMP * TEMP
              A(I, N) = CTEMP * A(I, N) - STEMP * TEMP
            end do
          end if
        end do
      else if (LSAME(DIRECT, 'B')) then
        do J = N - 1, 1, -1
          CTEMP = C(J)
          STEMP = S(J)
          if ((CTEMP /= ONE) .or. (STEMP /= ZERO)) then
            do I = 1, M
              TEMP = A(I, J)
              A(I, J) = STEMP * A(I, N) + CTEMP * TEMP
              A(I, N) = CTEMP * A(I, N) - STEMP * TEMP
            end do
          end if
        end do
      end if
    end if
  end if
  !
  return
  !
  !end of SLASR
  !
end
