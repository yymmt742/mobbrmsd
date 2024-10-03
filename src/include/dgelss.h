      pure subroutine DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, &
     &                        RANK, WORK, LWORK, INFO )
      use LA_CONSTANTS, only: wp=>dp
      integer,intent(in)      :: LDA, LDB, LWORK, M, N, NRHS
      integer,intent(out)     :: INFO, RANK
      real(wp),intent(in)     :: RCOND
      real(wp),intent(inout)  :: A( LDA, * ), B( LDB, * )
      real(wp),intent(out)    :: S( * ), WORK( * )
      end subroutine DGELSS
