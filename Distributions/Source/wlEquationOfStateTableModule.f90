MODULE wlEquationOfStateTableModule

  USE wlKindModule, ONLY: dp
  USE HDF5
  USE wlThermoStateModule
  USE wlDependentVariablesModule

  implicit none
  PRIVATE

  TYPE, PUBLIC :: EquationOfStateTableType
    !CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Name
    INTEGER                                      :: nVariables
    INTEGER, DIMENSION(3)                        :: nPoints
    TYPE(ThermoStateType)           :: TS 
    TYPE(DependentVariablesType)    :: DV
  END TYPE


  PUBLIC AllocateEquationOfStateTable
!  PUBLIC DeAllocateEquationOfStateTable

CONTAINS 

  SUBROUTINE AllocateEquationOfStateTable( EOSTable, nPoints, nVariables )

    TYPE(EquationOfStateTableType), INTENT(inout)    :: EOSTable
    INTEGER, INTENT(in)               :: nVariables
    INTEGER, DIMENSION(3), INTENT(in) :: nPoints
   
    EOSTable % nPoints(1:3) = nPoints(1:3)
    EOSTable % nVariables = nVariables

    CALL AllocateThermoState( EOSTable % TS, EOSTable % nPoints )
    CALL AllocateDependentVariables( EOSTable % DV, EOSTable % nPoints, &
                                     EOSTable % nVariables ) 

  END SUBROUTINE AllocateEquationOfStateTable

END MODULE wlEquationOfStateTableModule
