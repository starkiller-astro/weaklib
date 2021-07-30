PROGRAM wlOpacityDerivative

  USE wlKindModule, ONLY: dp
  USE HDF5
  
  USE wlIOModuleHDF
  USE wlEnergyGridModule
  USE wlOpacityTableModule
  USE wlOpacityFieldsModule
  USE wlOpacityTableIOModuleHDF
  USE wlExtNumericalModule, ONLY: epsilon
  USE wlInterpolationModule 
implicit none
    
   TYPE(OpacityTableType)  :: OpTD, OpacityTable
   INTEGER                 :: nOpacA = 4
   INTEGER                 :: nOpacB = 0
   INTEGER                 :: nMomB  = 0
   INTEGER                 :: nOpacC = 0
   INTEGER                 :: nMomC  = 0
   INTEGER                 :: nPointsE = 40   

   INTEGER                 :: i_r, i_e, j_rho, k_t, l_ye 
   REAL(dp), DIMENSION(1,4)  :: bufferD
   REAL(dp), DIMENSION(1)  :: x1, x2, x3, x4, temp
PRINT*, "Allocating OpacityTable"   

   CALL InitializeHDF( ) 
   CALL AllocateOpacityTable &
            ( OpTD, nOpacA, nOpacB, nMomB, &
              nOpacC, nMomC, nPointsE ) 
   CALL FinalizeHDF( )
 
   OpTD % thermEmAb % nOpacities = nOpacA

   OpTD % thermEmAb % nPoints(1) = nPointsE
   OpTD % thermEmAb % nPoints(2:4) = OpTD % nPointsTS

   OpTD % thermEmAb % Names = &
                                (/'Derivative E    ',&
                                  'Derivative rho  ',&
                                  'Derivative T    ',&
                                  'Derivative Ye   '/)
   OpTD % thermEmAb % Species = &
                                (/'Electron Neutrino', &
                                  'Electron Neutrino', &
                                  'Electron Neutrino', &
                                  'Electron Neutrino' /)
   OpTD % thermEmAb % Units = &
                                (/'N/A              ', &
                                  'N/A              ', &
                                  'N/A              ', &
                                  'N/A              '/)
   OpTD % thermEmAb % Offset = 1.0d-100
 
!-------------------------------------------------------------------------
WRITE (*,*) "Reading Opacity"
  CALL InitializeHDF( )  
  CALL ReadOpacityTableHDF( OpacityTable, 'OpacityTable.h5' )
  CALL FinalizeHDF( )

!-------------------------------------------------------------------------
!              Fill OpacityTable
!-------------------------------------------------------------------------
WRITE (*,*) "Fill OpTD table"
  ASSOCIATE( Deri => OpTD % thermEmAb % Absorptivity, &
            Table => OpacityTable % thermEmAb % Absorptivity(1) % Values, &
              EOS => OpacityTable % EOSTable % TS % States )
     DO l_ye = 1, OpacityTable % nPointsTS(3)
        x4(1) = EOS(3) % Values(l_ye)
       DO k_t = 1, OpacityTable % nPointsTS(2)
          x3(1) = EOS(2) % Values(k_t)
         DO j_rho = 1, OpacityTable % nPointsTS(1)
            x2(1) = EOS(1) % Values(j_rho)
           DO i_e = 1, OpacityTable % nPointsE
             x1(1) = OpacityTable % EnergyGrid % Values(i_e)
             CALL LogInterpolateDifferentiateSingleVariable&
             ( x1, x2, x3, x4, &
             OpacityTable % EnergyGrid % Values, &
             OpacityTable % EOSTable % TS % States(1) % Values, &
             OpacityTable % EOSTable % TS % States(2) % Values, &
             OpacityTable % EOSTable % TS % States(3) % Values, &
             (/1,1,1,0/), epsilon, Table, temp, bufferD, .FALSE. )
  
             DO i_r = 1, nOpacA            
               Deri( i_r ) % Values (i_e, j_rho, k_t, l_ye) = bufferD(1,i_r)
             END DO  !i_r

           END DO  !i_e
         END DO  !j_rho
       END DO  !k_t
     END DO  !l_ye
  END ASSOCIATE
 
WRITE (*,*) 'Describe OpTD '
  CALL DescribeOpacityTable( OpTD )

  CALL InitializeHDF( )
  CALL WriteOpacityTableHDF( OpTD, "OpacityDerivative.h5" )
  CALL FinalizeHDF( )

  WRITE (*,*) "HDF OpacityDerivative.h5 write successful"

!=============================================================

END PROGRAM wlOpacityDerivative
