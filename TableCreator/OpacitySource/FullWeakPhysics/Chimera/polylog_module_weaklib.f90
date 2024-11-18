MODULE polylog_module_weaklib

  USE wlKindModule, ONLY: dp

  implicit none

  integer, parameter :: nt=1000
  real(dp)   :: xa(0:nt), pl2a(0:nt), pl3a(0:nt) !Polylog table

end MODULE polylog_module_weaklib
