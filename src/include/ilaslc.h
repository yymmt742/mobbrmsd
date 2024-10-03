pure function ILASLC(M, N, A, LDA)
  integer, intent(in) :: M, N, LDA
  real, intent(in)    :: A(LDA, *)
  integer             :: ILASLC
end function ILASLC
