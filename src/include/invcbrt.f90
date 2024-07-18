subroutine invcbrt(x)
  implicit none
  real    :: x, fx, c, y
  integer :: iz
  x = 6.3 * 10
  fx = x / 3
  iz = INT(z"54a21d2a") - TRANSFER(x, iz) / 3
  !iz = INT(z"548c39cb") - TRANSFER(x, iz) / 3
  !iz = INT(z"54a21d2a") - TRANSFER(x, iz) / 3
  !print'(b64)', TRANSFER(ISHFT(TRANSFER(x, iz) * INT(z"55555556", INT64), -32), 0_INT32)
  !print'(b64)', INT(z"54A21D2A", INT32)
  !iz = INT(z"54A21D2A") - TRANSFER(x, iz) / 3
  !iz = INT(z"54A21D2A", INT32) - TRANSFER(ISHFT(TRANSFER(x, iz) * INT(z"55555556", INT64), -32), 0_INT32)
  print'(b64)', iz
  y = TRANSFER(iz, y)
  print'(b64.64)', iz
  print'(b64.64)', TRANSFER(y, iz)
  print *, 0, y, y - 1 / (x**(1 / 3.0))
  c = y * y * y * x
  !y = y * (1.5555555555 - c * (0.7777777777 - 0.2222222222 * c))
  y = y * (4.0 / 3 - y * y * y * fx)
  print *, 1, y, y**3, 1 / y**3, y - 1 / (x**(1 / 3.0))
  c = 1.0 - y * y * y * x
  y = y * (4.0 / 3 - y * y * y * fx)
  !y = y * (1.0000000000 - 0.3333333333 * c)
  print *, 2, y, y**3, 1 / y**3, y - 1 / (x**(1 / 3.0))
  !y = y * (1.5015480449 - 0.534850249 * y * y * y * x)
  !print *, y, y - 1 / (x**(1 / 3.0))
  !y = y * (1.333333985 - 0.33333333 * y * y * y * x)
  y = y * (4.0 / 3 - y * y * y * fx)
  print *, 3, y, y**3, 1 / y**3, y - 1 / (x**(1 / 3.0))
  y = y * (4.0 / 3 - y * y * y * fx)
  print *, 4, y, y**3, 1 / y**3, y - 1 / (x**(1 / 3.0))
  y = y * (4.0 / 3 - y * y * y * fx)
  print *, 5, y, y**3, 1 / y**3, y - 1 / (x**(1 / 3.0))
  end subroutine invcbrt(x)
