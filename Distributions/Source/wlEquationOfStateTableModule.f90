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
  PUBLIC DeAllocateEquationOfStateTable
  PUBLIC TableLimitFail

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

  SUBROUTINE DeAllocateEquationOfStateTable( EOSTable )

    TYPE(EquationOfStateTableType) :: EOSTable

    CALL DeAllocateThermoState( EOSTable % TS )
    CALL DeAllocateDependentVariables( EOSTable % DV )

  END SUBROUTINE DeAllocateEquationOfStateTable

  LOGICAL FUNCTION TableLimitFail( rho, t, ye, EOSTable )

    !LOGICAL                                    :: TableLimitFail
    REAL(dp), INTENT(in)                       :: rho, t, ye
    TYPE(EquationOfStateTableType), INTENT(in) :: EOSTable

      TableLimitFail = .false.
      IF ( rho < EOSTable % TS % minValues(1) ) TableLimitFail = .true.
      IF ( rho > EOSTable % TS % maxValues(1) ) TableLimitFail = .true.
      IF (   t < EOSTable % TS % minValues(2) ) TableLimitFail = .true.
      IF (   t > EOSTable % TS % maxValues(2) ) TableLimitFail = .true.
      IF (  ye < EOSTable % TS % minValues(3) ) TableLimitFail = .true.
      IF (  ye > EOSTable % TS % maxValues(3) ) TableLimitFail = .true.

  END FUNCTION TableLimitFail

END MODULE wlEquationOfStateTableModule
