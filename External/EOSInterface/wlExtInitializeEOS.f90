SUBROUTINE wlExtInitializeEOS( LSFilePath, LScompress )

!-- Initialize Equations of state and load required input data files

  CHARACTER(len=*), INTENT(in) :: LSFilePath
  CHARACTER(len=3), INTENT(in) :: LScompress

!-- Initialize BCK equation of state

  CALL eos0

!-- Initialize Lattimer-Swesty equation of state

  CALL loadmx( LSFilePath, LScompress )

END SUBROUTINE wlExtInitializeEOS
