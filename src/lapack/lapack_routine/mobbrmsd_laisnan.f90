  pure elemental function LA_ISNAN(x)
    real(RK), intent(in) :: x
    logical              :: LA_ISNAN
    LA_ISNAN = (x /= x)
  end function LA_ISNAN

