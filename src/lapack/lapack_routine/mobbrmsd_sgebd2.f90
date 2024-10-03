!> \brief \b SGEBD2 reduces a general matrix to bidiagonal form using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEBD2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgebd2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgebd2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgebd2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), D( * ), E( * ), TAUP( * ),
!      $                   TAUQ( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEBD2 reduces a real general m by n matrix A to upper or lower
!> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
!>
!> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows in the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns in the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the m by n general matrix to be reduced.
!>          On exit,
!>          if m >= n, the diagonal and the first superdiagonal are
!>            overwritten with the upper bidiagonal matrix B; the
!>            elements below the diagonal, with the array TAUQ, represent
!>            the orthogonal matrix Q as a product of elementary
!>            reflectors, and the elements above the first superdiagonal,
!>            with the array TAUP, represent the orthogonal matrix P as
!>            a product of elementary reflectors;
!>          if m < n, the diagonal and the first subdiagonal are
!>            overwritten with the lower bidiagonal matrix B; the
!>            elements below the first subdiagonal, with the array TAUQ,
!>            represent the orthogonal matrix Q as a product of
!>            elementary reflectors, and the elements above the diagonal,
!>            with the array TAUP, represent the orthogonal matrix P as
!>            a product of elementary reflectors.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is REAL array, dimension (min(M,N))
!>          The diagonal elements of the bidiagonal matrix B:
!>          D(i) = A(i,i).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is REAL array, dimension (min(M,N)-1)
!>          The off-diagonal elements of the bidiagonal matrix B:
!>          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
!>          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
!> \endverbatim
!>
!> \param[out] TAUQ
!> \verbatim
!>          TAUQ is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors which
!>          represent the orthogonal matrix Q. See Further Details.
!> \endverbatim
!>
!> \param[out] TAUP
!> \verbatim
!>          TAUP is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors which
!>          represent the orthogonal matrix P. See Further Details.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (max(M,N))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit.
!>          < 0: if INFO = -i, the i-th argument had an illegal value.
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
!> \date June 2017
!
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrices Q and P are represented as products of elementary
!>  reflectors:
!>
!>  If m >= n,
!>
!>     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
!>
!>  Each H(i) and G(i) has the form:
!>
!>     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!>
!>  where tauq and taup are real scalars, and v and u are real vectors;
!>  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
!>  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
!>  tauq is stored in TAUQ(i) and taup in TAUP(i).
!>
!>  If m < n,
!>
!>     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
!>
!>  Each H(i) and G(i) has the form:
!>
!>     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!>
!>  where tauq and taup are real scalars, and v and u are real vectors;
!>  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
!>  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
!>  tauq is stored in TAUQ(i) and taup in TAUP(i).
!>
!>  The contents of A on exit are illustrated by the following examples:
!>
!>  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!>
!>    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
!>    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
!>    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
!>    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
!>    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
!>    (  v1  v2  v3  v4  v5 )
!>
!>  where d and e denote diagonal and off-diagonal elements of B, vi
!>  denotes an element of the vector defining H(i), and ui an element of
!>  the vector defining G(i).
!> \endverbatim
!>
!  =====================================================================
pure subroutine SGEBD2(M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO)
  implicit none
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
  integer, intent(in)  :: LDA, M, N
  integer, intent(out) :: INFO
!     ..
!     .. Array Arguments ..
  real(RK), intent(inout)  :: A(LDA, *)
  real(RK), intent(out)    :: D(*), E(*), TAUP(*), TAUQ(*), WORK(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
  integer :: I
!
!     .. Parameters ..
! real(RK), parameter      :: ZERO = 0.0E+0, ONE = 1.0E+0
!     ..
! interface
! .. External Subroutines ..
!   include 'slarf.h'
!   include 'slarfg.h'
!   include 'xerbla.h'
! end interface
!..
!..intrinsic Functions..
  intrinsic :: MAX, MIN
!..
!..Executable Statements..
!
! Test the input parameters
!
  INFO = 0
  if (M < 0) then
    INFO = -1
  else if (N < 0) then
    INFO = -2
  else if (LDA < MAX(1, M)) then
    INFO = -4
  end if
  if (INFO < 0) then
!   call XERBLA('SGEBD2', -INFO)
    return
  end if
!
  if (M >= N) then
!
! Reduce to upper bidiagonal form
!
    do I = 1, N
!
! Generate elementary reflector H(i) to annihilate A(i + 1:m, i)
!
      call SLARFG(M - I + 1, A(I, I), A(MIN(I + 1, M), I), 1, TAUQ(I))
      D(I) = A(I, I)
      A(I, I) = ONE
!
! Apply H(i) to A(i:m, i + 1:n) from the left
!
      if (I < N) call SLARF('Left', M - I + 1, N - I, A(I, I), 1, TAUQ(I), A(I, I + 1), LDA, WORK)
      A(I, I) = D(I)
!
      if (I < N) then
!
! Generate elementary reflector G(i) to annihilate A(i, i + 2:n)
!
        call SLARFG(N - I, A(I, I + 1), A(I, MIN(I + 2, N)), LDA, TAUP(I))
        E(I) = A(I, I + 1)
        A(I, I + 1) = ONE
!
! Apply G(i) to A(i + 1:m, i + 1:n) from the right
!
        call SLARF('Right', M - I, N - I, A(I, I + 1), LDA, TAUP(I), A(I + 1, I + 1), LDA, WORK)
        A(I, I + 1) = E(I)
      else
        TAUP(I) = ZERO
      end if
    end do
  else
!
! Reduce to lower bidiagonal form
!
    do I = 1, M
!
! Generate elementary reflector G(i) to annihilate A(i, i + 1:n)
!
      call SLARFG(N - I + 1, A(I, I), A(I, MIN(I + 1, N)), LDA, TAUP(I))
      D(I) = A(I, I)
      A(I, I) = ONE
!
! Apply G(i) to A(i + 1:m, i:n) from the right
!
      if (I < M) call SLARF('Right', M - I, N - I + 1, A(I, I), LDA, TAUP(I), A(I + 1, I), LDA, WORK)
      A(I, I) = D(I)
!
      if (I < M) then
!
! Generate elementary reflector H(i) to annihilate A(i + 2:m, i)
!
        call SLARFG(M - I, A(I + 1, I), A(MIN(I + 2, M), I), 1, TAUQ(I))
        E(I) = A(I + 1, I)
        A(I + 1, I) = ONE
!
! Apply H(i) to A(i + 1:m, i + 1:n) from the left
!
        call SLARF('Left', M - I, N - I, A(I + 1, I), 1, TAUQ(I), A(I + 1, I + 1), LDA, WORK)
        A(I + 1, I) = E(I)
      else
        TAUQ(I) = ZERO
      end if
    end do
  end if
  return
!
! end of SGEBD2
!
end
