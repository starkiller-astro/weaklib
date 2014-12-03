PROGRAM wlStateModuleTest

  USE wlStateModule
!  USE wlKindModule, ONLY: dp

  implicit none


INTEGER, DIMENSION(3) :: npts
INTEGER :: nvar
INTEGER :: j   
TYPE(StateType) :: State1

  npts = (/2,2,2/)


  Call AllocateState( State1, npts )

  State1%nValues(1:3) = npts(1:3)
  State1%Names(1:3) = (/'Density                         ',&
                        'Temperature                     ',&
                        'Electron Fraction               '/)
    
  State1%minValues(1:3) =  (/1.0d06,0.1d00,1.0d-02/)
  State1%maxValues(1:3) =  (/1.0d15,5.0d01,6.1d-01/)

  DO nvar = 1,3
    DO j = 1,State1%nValues(nvar)
      State1%States(nvar)%Values(1:j) = State1%minValues(1:j) + (10.d0**((DBLE(j-1)/DBLE(State1%nValues(1:j)))*LOG10(State1%maxValues(1:j)-State1%minValues(1:j)))) 
    END DO
  END DO

  WRITE(*,*) State1%States(1)%Values(1)
  WRITE(*,*) State1%States(1)%Values(2)
  WRITE(*,*) State1%States(1)%Values(3)
  WRITE(*,*) State1%States(2)%Values(1)
  WRITE(*,*) State1%States(2)%Values(2)
  WRITE(*,*) State1%States(2)%Values(3)
  WRITE(*,*) State1%States(3)%Values(1)
  WRITE(*,*) State1%States(3)%Values(2)
  WRITE(*,*) State1%States(3)%Values(3)


!  Call DeAllocateState( State1, npts )

  !TYPE(State) :: EOS
  !EOS = State( (/2,2,2/) )

END PROGRAM wlStateModuleTest
