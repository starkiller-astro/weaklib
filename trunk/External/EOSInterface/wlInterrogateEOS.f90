PROGRAM wlInterrogateEOS

!-- User interactive interrogation of external equations of state

USE wlKindModule, ONLY: dp
USE wlExtEOSWrapperModule, ONLY: wlGetFullEOS

  CHARACTER(len=1)   :: EOSFlag     ! nuclear eos selection flag
  CHARACTER(len=3)   :: LScompress
  CHARACTER(len=128) :: LSFilePath

  REAL(dp)           :: Density     ! Density [g/cm3]
  REAL(dp)           :: Temperature ! Temperature [K]
  REAL(dp)           :: Ye          ! Electron fraction

  LOGICAL            :: fail        ! did EoS fail to converge

  REAL(dp)           :: press       ! pressure
  REAL(dp)           :: energ       ! internal energy
  REAL(dp)           :: entrop      ! entropy [kb/baryon]
  REAL(dp)           :: chem_n      ! free neutron cemical potential
  REAL(dp)           :: chem_p      ! free proton chemical potential
  REAL(dp)           :: chem_e      ! electron chemical potential
  REAL(dp)           :: xn_neut     ! free neutron fraction
  REAL(dp)           :: xn_prot     ! free proton fraction
  REAL(dp)           :: xn_alpha    ! alpha fraction
  REAL(dp)           :: xn_heavy    ! heavy fraction
  REAL(dp)           :: a_heavy     ! A for mean heavy nucleus
  REAL(dp)           :: z_heavy     ! Z for mean heavy nucleus
  REAL(dp)           :: be_heavy    ! Binding energy for mean heavy nucleus
  REAL(dp)           :: thermalenergy ! Internal energy without rest mass and constant offsets
  REAL(dp)           :: gamma1 

  INTEGER            :: i

  CHARACTER(len=32), DIMENSION(14) :: Names
  CHARACTER(len=32), DIMENSION(14) :: Units

  REAL(dp), DIMENSION(14)          :: Values

 99  FORMAT(a32,'= ', es14.6, ' ', a32)

!-- Initialize EoSs and labels

  LScompress = '220'
  LSFilePath = '../../LS/Data'
  CALL  wlExtInitializeEOS( LSFilePath, LScompress )

  Names(1:15) = (/'Pressure                        ', &
                  'Entropy Per Baryon              ', &
                  'Internal Energy Density         ', &
                  'Neutron Chemical Potential      ', &
                  'Electron Chemical Potential     ', &
                  'Proton Chemical Potential       ', &
                  'Neutron Mass Fraction           ', &
                  'Proton Mass Fraction            ', &
                  'Alpha Mass Fraction             ', &
                  'Heavy Mass Fraction             ', &
                  'Heavy Mass Number               ', &
                  'Heavy Charge Number             ', &
                  'Heavy Binding Energy            ', &
                  'Thermal Energy                  ', &
                  'Gamma1                          '/)

  Units(1:15) = (/'Dynes per cm^2                  ', &
                  'k_b per baryon                  ', &
                  'erg per gram                    ', &
                  'MeV                             ', &
                  'MeV                             ', &
                  'MeV                             ', &
                  '                                ', &
                  '                                ', &
                  '                                ', &
                  '                                ', &
                  '                                ', &
                  '                                ', &
                  'MeV                             ', &
                  'MeV                             ', &
                  '                                '/)

!-- Take user input until interrupted

  DO
    WRITE(*,*) 'Enter input parameters for EOS (EOS selection["L" or "B", density [g/cm^3], temperature[K], Ye) '
    WRITE(*,*) 'example: "L",1.5e13, 3.4e10, 0.4'
    READ(*,*) EOSFlag, Density, Temperature, Ye
    CALL wlGetFullEOS( Density, Temperature, Ye, EOSFlag, fail, press, energ, &
                       entrop, chem_n, chem_p, chem_e, xn_neut, xn_prot,      &
                       xn_alpha, xn_heavy, a_heavy, z_heavy, be_heavy, thermalenergy, gamma1 )

    Values(1)  = press   
    Values(3)  = energ   
    Values(2)  = entrop  
    Values(4)  = chem_n  
    Values(6)  = chem_p  
    Values(5)  = chem_e  
    Values(7)  = xn_neut 
    Values(8)  = xn_prot 
    Values(9)  = xn_alpha
    Values(10) = xn_heavy
    Values(11) = a_heavy 
    Values(12) = z_heavy 
    Values(13) = be_heavy
    Values(14) = thermalenergy
    Values(15) = gamma1

    WRITE(*,*) "EOS ", EOSFlag, " returns dependent variables:"

    DO i = 1,15
      WRITE(*,99) Names(i), Values(i), ADJUSTL( Units(i) )
    END DO

  END DO

END PROGRAM wlInterrogateEOS
