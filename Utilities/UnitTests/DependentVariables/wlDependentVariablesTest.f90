PROGRAM wlDependentVariablesTest
 
  USE wlDependentVariablesModule
 
  implicit none

  INTEGER :: i
  TYPE(DependentVariablesType) :: DV
 
  CALL AllocateDependentVariables( DV, nVariables = 5, nValues = (/10,10,10/) )

  DO i = 1, SIZE( DV % Variables )
    WRITE (*,*) SHAPE( DV % Variables(i) % Values )
  END DO

  CALL DeAllocateDependentVariables( DV )
 
END PROGRAM wlDependentVariablesTest
