MODULE polylog_module_weaklib

  USE kind_module, ONLY: double

  implicit none

  integer, parameter :: nt=1000
  real(double)   :: xa(0:nt), pl2a(0:nt), pl3a(0:nt) !Polylog table

end MODULE polylog_module_weaklib
