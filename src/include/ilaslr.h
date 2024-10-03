pure function ILASLR(M, N, A, LDA)
  integer, intent(in) :: M, N, LDA
  real, intent(in)    :: A(LDA, *)
  integer             :: ILASLR
end function ILASLR
