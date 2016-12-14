PROGRAM wlProfileInterpolationTest

  USE wlKindModule, ONLY: dp 
  USE wlGridModule
  USE HDF5
  USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS
  USE wlExtPhysicalConstantsModule
  USE wlEquationOfStateTableModule
  USE wlInterpolationModule
  USE wlIOModuleHDF
  USE wlEOSIOModuleHDF
  USE sfho_frdm_composition_module
  !USE sfhx_frdm_composition_module

  implicit none

  INTEGER  :: i, j, k, TestUnit1, TestUnit2, TestUnit3, FileUnit 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rho
  REAL(dp), DIMENSION(:), ALLOCATABLE :: T
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Ye 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: heavymassfrac, heavyA, heavyZ, heavyBE
  INTEGER :: NumGoodPoints
  TYPE(EquationOfStateTableType) :: EOSTable

  REAL(dp), DIMENSION(:), ALLOCATABLE :: Interpolant
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: DirectCall 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Interpolants 
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Derivative
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Derivatives 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: press 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: entrop
  REAL(dp), DIMENSION(:), ALLOCATABLE :: energ
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_e
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_p
  REAL(dp), DIMENSION(:), ALLOCATABLE :: chem_n 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn_prot 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn_neut 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn_alpha 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xn_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: z_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: a_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: be_heavy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: thermenergy 
  REAL(dp), DIMENSION(:), ALLOCATABLE :: gammaone
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rand
  REAL(dp) :: Yemin, Yemax, logrho, logTmin, logTmax, logrhomin, logrhomax, epsilon
  INTEGER, PARAMETER :: NumPoints = 292 
  REAL(dp), DIMENSION(15) :: L1norm, Maxnorm, L1norm2, Maxnorm2

  REAL(dp)                              :: tb,yeb,nb

!.....output:
!.....naz(inuc): number density of nucleus 'inuc' [fm^-3]
!.....nn: number density of neutrons [fm^-3]
!.....np: number density of protons [fm^-3]
!.....xaz(inuc): mass fraction of nucleus 'inuc' []
!.....xn: mass fraction of neutrons []
!.....xp: mass fraction of protons []
!.....sflag: success flag
!.....sflag=0: input is not inside the table and the composition cannot
!.....         be calculated
!.....sflag=1: succesful call of the routine

   REAL(dp), DIMENSION(kmax)              :: naz, xaz
   REAL(dp) :: xn, xp, nn, np
   INTEGER :: sflag

  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath

  LOGICAL            :: fail        ! did EoS fail to converge

  99 FORMAT ("rho=", es12.5, " T=", es12.5, " Ye=" , es12.5, " Int=", es12.5, &
             " DC=", es12.5, " Diff=", es12.5 ) 
 4018 FORMAT (i4,3(es12.5))
  epsilon = 1.d-100


  FileUnit = 50 
  OPEN( unit = FileUnit, FILE="Output0ms.d")
  OPEN( newunit = TestUnit1, FILE="InterpolateAllTest.d")
  !OPEN( newunit = TestUnit2, FILE="EOSComparisonTestSHFx.txt")
  OPEN( newunit = TestUnit2, FILE="EOSComparisonTestSFHoBEtest.d")
  !OPEN( newunit = TestUnit2, FILE="EOSComparisonTestSTOS.d")
  !OPEN( newunit = TestUnit2, FILE="EOSComparisonTestCLS.d")
  OPEN( newunit = TestUnit3, FILE="HeavyComparisonTestSFHoBEtest.d")

  LScompress = '220'
  LSFilePath = '../../../LS/Data'
  EOSFlag = "B"

  CALL InitializeHDF( )

  !CALL ReadEquationOfStateTableHDF( EOSTable, "wl-EOS-LS220-20-40-100.h5" )
  !CALL ReadEquationOfStateTableHDF( EOSTable, "EquationOfStateTable.h5" )
  !CALL ReadEquationOfStateTableHDF( EOSTable, "WeakLibCLS.h5" )
  !CALL ReadEquationOfStateTableHDF( EOSTable, "wl-EOS-SFHx-25-25-100.h5" )
  !CALL ReadEquationOfStateTableHDF( EOSTable, "wl-EOS-SFHo-25-25-100.h5" )
  CALL ReadEquationOfStateTableHDF( EOSTable, "SFHoTEST.h5" )
  !CALL ReadEquationOfStateTableHDF( EOSTable, "wl-EOS-STOS-25-25-100.h5" )

  ALLOCATE( rho( NumPoints ), T( NumPoints ), Ye( NumPoints ),  &
            DirectCall( NumPoints, 15), press( NumPoints ), &
            entrop( NumPoints ), energ( NumPoints ), chem_e( NumPoints ),                & 
            chem_p( NumPoints ), chem_n( NumPoints ), xn_prot( NumPoints ),              &
            xn_neut( NumPoints ), xn_alpha( NumPoints ), xn_heavy( NumPoints ),          &
            z_heavy( NumPoints ), a_heavy( NumPoints ), be_heavy( NumPoints ),           &
            thermenergy( NumPoints ), gammaone( NumPoints ), heavymassfrac( NumPoints ), &
            Interpolants( NumPoints, EOSTable % DV % nVariables  ), heavyZ( NumPoints ), &
            heavyA( NumPoints), heavyBE( NumPoints) )                   


  DO i = 1, NumPoints
    READ(FileUnit, '(3(es11.3,x))') rho(i), T(i), Ye(i)
  END DO

!  CALL MakeLogGrid(1.05d08, 9.95d11,NumPoints,rho)
!  T(1:NumPoints) = 1.0d11
!  Ye(1:NumPoints) = 4.0d-01

  CALL wlExtInitializeEOS( LSFilePath, LScompress )

  DO i = 1, SIZE( rho )

    CALL wlGetFullEOS( rho(i), T(i), Ye(i), EOSFlag, fail, press(i), energ(i), &
                       entrop(i), chem_n(i), chem_p(i), chem_e(i), xn_neut(i), &
                       xn_prot(i), xn_alpha(i), xn_heavy(i), a_heavy(i),       &
                       z_heavy(i), be_heavy(i), thermenergy(i), gammaone(i), 1, 1, 1 ) 

    DirectCall(i,1) = press(i)
    DirectCall(i,2) = entrop(i)
    DirectCall(i,3) = energ(i)
    DirectCall(i,4) = chem_e(i)
    DirectCall(i,5) = chem_p(i)
    DirectCall(i,6) = chem_n(i)
    DirectCall(i,7) = xn_neut(i)
    DirectCall(i,8) = xn_prot(i)
    DirectCall(i,9) = xn_alpha(i)
    DirectCall(i,10) = xn_heavy(i)
    DirectCall(i,11) = a_heavy(i)
    DirectCall(i,12) = z_heavy(i)
    DirectCall(i,13) = be_heavy(i)
    DirectCall(i,14) = thermenergy(i)
    DirectCall(i,15) = gammaone(i)

  END DO
  
  NumGoodPoints = 0
  
  CALL LogInterpolateAllVariables( rho, T, Ye, EOSTable % TS % LogInterp, &
         EOSTable % TS, EOSTable % DV, Interpolants ) 
    DO i = 1, SIZE(rho)
      DO j = 1, EOSTable % DV % nVariables 
        L1norm(j) = ABS(( DirectCall(i,j) - Interpolants(i,j) )/ DirectCall(i,j))
        WRITE(*,*) i, rho(i), DirectCall(i,j), Interpolants(i,j) 
      END DO

      WRITE (TestUnit1,'( 15(es14.7,x) )' ) L1norm(1), L1norm(2),   &
                                  L1norm(3), L1norm(4),   &
                                  L1norm(5), L1norm(6),   &
                                  L1norm(7), L1norm(8),   &
                                  L1norm(9), L1norm(10),  &
                                  L1norm(11), L1norm(12), &
                                  L1norm(13), L1norm(14), &
                                  L1norm(15)
      !z(8077,2)*naz(8077)WRITE (TestUnit2, '(3(es14.7,x))' ) rho(i), DirectCall(i,6), Interpolants(i,6) 
      WRITE (TestUnit2, '(18(es14.7,x))' ) rho(i), T(i), Ye(i), &
                                          Interpolants(i,1), Interpolants(i,2), Interpolants(i,3), & 
                                          Interpolants(i,4), Interpolants(i,5), Interpolants(i,6), & 
                                          Interpolants(i,7), Interpolants(i,8), Interpolants(i,9), & 
                                          Interpolants(i,10), Interpolants(i,11), Interpolants(i,12), & 
                                          Interpolants(i,13), Interpolants(i,14), Interpolants(i,15)
    END DO

    CALL compdata_readin

    DO i = 1, SIZE(rho)
      tb=T(i) * kmev 
      yeb=Ye(i)
      nb=rho(i) * kfm
      CALL sub_dist_interpol(tb,yeb,nb,xaz,xn,xp,naz,nn,np,sflag)
      heavymassfrac(i) = SUM( az(:,1)*naz/nb) - az(8077,1)*naz(8077)/nb &
                              - az(8076,1)*naz(8076)/nb - az(8075,1)*naz(8075)/nb &
                              - az(8074,1)*naz(8074)/nb
      IF ( heavymassfrac(i) == 0.0d0 ) THEN
        heavyZ(i) = 1.8d0
      ELSE
      heavyZ(i) = ( SUM( az(:,2)*naz) - az(8077,2)*naz(8077) - az(8076,2)*naz(8076) &
                    - az(8075,2)*naz(8075) - az(8074,2)*naz(8074) )/( SUM( naz ) - naz(8077) - naz(8076) &
                    - naz(8075) - naz(8074))
      END IF
      IF ( heavymassfrac(i) == 0.0d0 ) THEN 
        heavyA(i) = 4.0d0
      ELSE
      heavyA(i) = ( SUM( az(:,1)*naz)- az(8077,1)*naz(8077) - az(8076,1)*naz(8076) &
                    - az(8075,1)*naz(8075) - az(8074,1)*naz(8074))/( SUM( naz ) - naz(8077) - naz(8076) &
                    - naz(8075) - naz(8074) )
      END IF
      IF ( heavymassfrac(i) == 0.0d0 ) THEN 
        heavyBE(i) = 0.0d0
      ELSE
        heavyBE(i) = ( SUM( -bind(:) * xaz(:) / az(:,1) )       & !+ bind(iHe4) * xaz(iHe4) / 4.d0   &     
             + bind(8077) * xaz(8077) / az(8077,1)    &
             + bind(8076) * xaz(8076) / az(8076,1)    &
             + bind(8075) * xaz(8075) / az(8075,1)    &
             + bind(8074) * xaz(8074) / az(8074,1) )/ heavymassfrac(i)
WRITE(*,*) 'heavyBE=', heavyBE(i)
      END IF
      WRITE (TestUnit3, '(14(es14.7,x))' ) rho(i), T(i), Ye(i), DirectCall(i,10), Interpolants(i,10), heavymassfrac(i), &
              DirectCall(i,11), Interpolants(i,11), heavyZ(i), DirectCall(i,12), Interpolants(i,12), heavyA(i), & 
              DirectCall(i,13), heavyBE(i)  
    END DO



  CLOSE(TestUnit1)
  CLOSE(TestUnit2)
  CLOSE(TestUnit3)
  CLOSE(FileUnit)

  CALL DeAllocateEquationOfStateTable( EOSTable )

  CALL FinalizeHDF( )

END PROGRAM wlProfileInterpolationTest
