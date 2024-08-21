MODULE wlGammaSoundSpeed
	
	USE wlKindModule, ONLY: dp
	USE wlExtNumericalModule, ONLY: zero, half, one, pi
	USE wlExtPhysicalConstantsModule, ONLY: kmev, rmu, kmev_inv, ergmev, me, cvel, cm3fm3
	
	IMPLICIT NONE
	PRIVATE
	
	PUBLIC :: derivatives_production
	
	CONTAINS
	
	subroutine derivatives_production(igamma, nrho, ntemp, nye, logrho, logtemp, ye, &
		eos_table, energy_shift, cs2, gamma)
		
		integer, intent(in) :: igamma
		integer, intent(in) :: nrho, ntemp, nye
		real(dp), intent(in) :: eos_table(nrho,ntemp,nye,3)
		real(dp), intent(in) :: energy_shift
		real(dp), intent(in) :: logrho(nrho), logtemp(ntemp), ye(nye)
		real(dp), allocatable, intent(out) :: cs2(:,:,:)
		real(dp), allocatable, intent(out) :: gamma(:,:,:)
		
		integer i,j,k
		real(dp) :: dx,x1,x2,f1,f2,z,zz,h,cp,cv,beta_v,gamma_ad,kappa_t
		
		real(dp), allocatable :: dedT(:,:,:)
		real(dp), allocatable :: dsdlnT(:,:,:)
		real(dp), allocatable :: dsdlnrho(:,:,:)
		real(dp), allocatable :: dpdrhoT(:,:,:)
		real(dp), allocatable :: dedrhoT(:,:,:)
		real(dp), allocatable :: dpdT(:,:,:)
		real(dp), allocatable :: dpdrho(:,:,:)
		real(dp), allocatable :: dpde(:,:,:)
		integer :: ipress = 1
		integer :: ientropy = 2
		integer :: ienergy = 3
		
		! igamma: There are multiple ways to compute gamma1, via the entropy
		!         turns out to be pretty good.
			! 1 -> Gamma1 = &
			!      dlnP/dlnrho|T,Y_e - ds/dlnrho|T,Y_e * (dlnP/dlnT / ds/dlnT)_rho,Y_e
			! 2 -> Gamma1 = P/rho * ( dpdrho|e,Y_e + P/rho**2 dpde|rho,Y_e
		! 3 -> Gamma1 = dlnP/dlnrho|T,Y_e + T * (dpdT|rho,Y_e)**2 / &
		!      (P*rho * dedT|rho,Y_e)
		
		!###########################################################################  
		! dedT, e is erg in log10, T is in MeV in log 10
		allocate(dedT(nrho,ntemp,nye))
		dedT = 0.0d0
		
		do k=1,nye
			do i=1,nrho
				do j=2,ntemp-1
					x1 = logtemp(j-1)
					f1 = eos_table(i,j-1,k,ienergy)
					x2 = logtemp(j+1)
					f2 = eos_table(i,j+1,k,ienergy)
					dedT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j)  &
					* 10.0d0**eos_table(i,j,k,ienergy)
				enddo
				
				! boundaries: one-sided derivative
				x2 = logtemp(2)
				x1 = logtemp(1)
				f2 = eos_table(i,2,k,ienergy)
				f1 = eos_table(i,1,k,ienergy)
				dedT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1)  &
                * 10.0d0**eos_table(i,1,k,ienergy)
				x2 = logtemp(ntemp)
				x1 = logtemp(ntemp-1)
				f2 = eos_table(i,ntemp,k,ienergy)
				f1 = eos_table(i,ntemp-1,k,ienergy)
				dedT(i,ntemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(ntemp) &
				* 10.0d0**eos_table(i,ntemp,k,ienergy) 
			enddo
		enddo
		
		!########################################################################### 
		
		!dpdT, p is in dyn/cm^2 in log10, T is in MeV in log 10
		!dsdlnT, s is in k_B/baryon not in log 10, T is in meV in log 10
		allocate(dpdT(nrho,ntemp,nye))
		allocate(dsdlnT(nrho,ntemp,nye))
		dpdT = 0.0d0
		dsdlnT = 0.0d0
		
		do k=1,nye
			do i=1,nrho
				do j=2,ntemp-1
					x1 = logtemp(j-1)
					f1 = eos_table(i,j-1,k,ipress)
					x2 = logtemp(j+1)
					f2 = eos_table(i,j+1,k,ipress)
					dpdT(i,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(j)  &
					* 10.0d0**eos_table(i,j,k,ipress)
					
					x1 = logtemp(j-1)
					f1 = log10(eos_table(i,j-1,k,ientropy))
					x2 = logtemp(j+1)
					f2 = log10(eos_table(i,j+1,k,ientropy))
					dsdlnT(i,j,k) = (f2-f1)/(x2-x1) * eos_table(i,j,k,ientropy)
					
				enddo
				
				! boundaries: one-sided derivative
				x1 = logtemp(1)
				f1 = eos_table(i,1,k,ipress)
				x2 = logtemp(2)
				f2 = eos_table(i,2,k,ipress)
				dpdT(i,1,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1) * &
				10.0d0**eos_table(i,1,k,ipress)
				
				x1 = logtemp(1)
				f1 = log10(eos_table(i,1,k,ientropy))
				x2 = logtemp(2)
				f2 = log10(eos_table(i,2,k,ientropy))
				dsdlnT(i,1,k) = (f2-f1)/(x2-x1) * eos_table(i,1,k,ientropy)
				
				
				x1 = logtemp(ntemp-1)
				f1 = eos_table(i,ntemp-1,k,ipress)
				x2 = logtemp(ntemp)
				f2 = eos_table(i,ntemp,k,ipress)
				dpdT(i,ntemp,k) = (f2-f1)/(x2-x1) / 10.0d0**logtemp(1) * &
				10.0d0**eos_table(i,1,k,ipress)
				
				x1 = logtemp(ntemp-1)
				f1 = log10(eos_table(i,ntemp-1,k,ientropy))
				x2 = logtemp(ntemp)
				f2 = log10(eos_table(i,ntemp,k,ientropy))
				dsdlnT(i,ntemp,k) = (f2-f1)/(x2-x1) & 
				* eos_table(i,ntemp,k,ientropy)
				
			enddo
		enddo
		
		!########################################################################### 
		! dp/drho|T, p is in dyn/cm^2 in log 10, T is in MeV in log 10
			! ds/dlnrho|T, s is in k_b/baryon not in log 10, T is in MeV in log 10
			! de/drho|T, e is in erg in log 10, T is in MeV in log 10
			! dp/drho|e = dp/drho|T + dp/dT * (-de/drho|T) / de/dT
		! dp/de|rho = dp/dT / de/dT
		
		allocate(dpdrho(nrho,ntemp,nye))
		allocate(dpde(nrho,ntemp,nye))
		allocate(dsdlnrho(nrho,ntemp,nye))
		allocate(dpdrhoT(nrho,ntemp,nye))
		allocate(dedrhoT(nrho,ntemp,nye))
		dpdrho   = 0.0d0
		dpde     = 0.0d0
		dsdlnrho = 0.0d0
		dpdrhoT  = 0.0d0
		dedrhoT  = 0.0d0
		
		do k=1,nye
			do j=1,ntemp
				do i=2,nrho-1
					x1 = logrho(i-1)
					x2 = logrho(i+1)
					f1 = eos_table(i-1,j,k,ipress)
					f2 = eos_table(i+1,j,k,ipress)
					dpdrhoT(i,j,k) = (f2-f1)/(x2-x1) /  10.0d0**logrho(i) &
					* 10.0d0**eos_table(i,j,k,ipress)
					
					x1 = logrho(i-1)
					x2 = logrho(i+1)
					f1 = eos_table(i-1,j,k,ienergy)
					f2 = eos_table(i+1,j,k,ienergy)
					dedrhoT(i,j,k) = (f2-f1)/(x2-x1) /  10.0d0**logrho(i) &
					* 10.0d0**eos_table(i,j,k,ienergy)
					
					x1 = logrho(i-1)
					x2 = logrho(i+1)
					f1 = log10(eos_table(i-1,j,k,ientropy))
					f2 = log10(eos_table(i+1,j,k,ientropy))
					dsdlnrho(i,j,k) = (f2-f1)/(x2-x1) * eos_table(i,j,k,ientropy)
				enddo
				
				! boundaries: one-sided derivative
				x1 = logrho(1)
				x2 = logrho(2)
				f1 = eos_table(1,j,k,ipress)
				f2 = eos_table(2,j,k,ipress)
				dpdrhoT(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
				* 10.0d0**eos_table(1,j,k,ipress)
				
				f1 = eos_table(1,j,k,ienergy)
				f2 = eos_table(2,j,k,ienergy)
				dedrhoT(1,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(1) &
				* 10.0d0**eos_table(1,j,k,ienergy)
				
				x1 = logrho(nrho-1)
				x2 = logrho(nrho)
				f1 = eos_table(nrho-1,j,k,ipress)
				f2 = eos_table(nrho,j,k,ipress)
				dpdrhoT(nrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nrho) &
				* 10.0d0**eos_table(nrho,j,k,ipress)
				
				f1 = eos_table(nrho-1,j,k,ienergy)
				f2 = eos_table(nrho,j,k,ienergy)
				dedrhoT(nrho,j,k) = (f2-f1)/(x2-x1) / 10.0d0**logrho(nrho) &
				* 10.0d0**eos_table(nrho,j,k,ienergy)
				
				x1 = logrho(1)
				x2 = logrho(2)
				f1 = log10(eos_table(1,j,k,ientropy))
				f2 = log10(eos_table(2,j,k,ientropy))
				dsdlnrho(1,j,k) = (f2-f1)/(x2-x1) * eos_table(1,j,k,ientropy)
				
				x1 = logrho(nrho-1)
				x2 = logrho(nrho)
				f1 = log10(eos_table(nrho-1,j,k,ientropy))
				f2 = log10(eos_table(nrho,j,k,ientropy))
				dsdlnrho(nrho,j,k) = (f2-f1)/(x2-x1) * eos_table(nrho,j,k,ientropy)
				
			enddo
		enddo
		
		do k=1,nye
			do j=1,ntemp
				do i=1,nrho
					dpdrho(i,j,k) = dpdrhoT(i,j,k) + dpdT(i,j,k) * &
					(-dedrhoT(i,j,k))/dedT(i,j,k)
					
					dpde(i,j,k) = dpdT(i,j,k) / dedT(i,j,k)
				enddo
			enddo
		enddo
		
		!########################################################################### 
		allocate(gamma(nrho,ntemp,nye))
		allocate(cs2(nrho,ntemp,nye))
		cs2 = 0.0d0
		gamma = 0.0d0
		
		do k=1,nye
			do j=1,ntemp
				do i=1,nrho
					z = 10.0d0**logrho(i)/10.0d0**eos_table(i,j,k,ipress)
					zz = 10.0d0**logtemp(j)/10.0d0**eos_table(i,j,k,ipress)
					
					if(igamma .eq. 1) then
						gamma(i,j,k) = z * dpdrhoT(i,j,k) & 
						- dsdlnrho(i,j,k) * zz * dpdT(i,j,k) / dsdlnT(i,j,k)
					endif
					
					if(igamma .eq. 2) then
						gamma(i,j,k) =  dpdrho(i,j,k) + 10.0d0**eos_table(i,j,k,ipress) / &
						(10.0d0**logrho(i))**2 * dpde(i,j,k)
						gamma(i,j,k) = z * gamma(i,j,k)
					endif
					
					if(igamma .eq. 3) then
						gamma(i,j,k) = 10.0d0**logrho(i)/10.0d0**eos_table(i,j,k,ipress) &
						* dpdrhoT(i,j,k) &
						+ 10.0d0**logtemp(j)*dpdT(i,j,k)**2 / &
						(10.0d0**eos_table(i,j,k,ipress)*10.0d0**logrho(i) * dedT(i,j,k))
					endif
					
					cs2(i,j,k) = gamma(i,j,k) / z
					cs2(i,j,k) = gamma(i,j,k) * 10.0d0**eos_table(i,j,k,ipress) / &
					(10.0d0**eos_table(i,j,k,ipress) + (10.0d0**eos_table(i,j,k,ienergy) - energy_shift - &
					0.511d0 / rmu * ergmev * ye(k) + cvel**2) * 10.0d0**logrho(i) ) * cvel**2
					
					if(igamma .eq. 5) then
						h = 10.0d0**eos_table(i,j,k,ipress) + (10.0d0**eos_table(i,j,k,ienergy) - energy_shift - &
							0.511d0 / rmu * ergmev * ye(k) + cvel**2) * 10.0d0**logrho(i)
						
						kappa_t = 10.0d0**logrho(i) * dpdrho(i,j,k) * (cm3fm3/ergmev)
						kappa_t = 1.0d0 / kappa_t
						beta_v = dpdT(i,j,k) * (cm3fm3/ergmev) / kmev
						cv = 1.0d0 / (10.0d0**logrho(i)/rmu*cm3fm3) * dsdlnT(i,j,k) / LOG(10.0d0) 
						cp = cv + 10.0d0**logtemp(j)*kmev/(10.0d0**logrho(i)/rmu*cm3fm3) * beta_v * kappa_t * beta_v
						gamma_ad = cp/cv
													
						cs2(i,j,k) = gamma_ad / (h * (cm3fm3/ergmev) * kappa_t )
						gamma(i,j,k) = gamma_ad / 10.0d0**eos_table(i,j,k,ipress)*ergmev/cm3fm3 / kappa_t
						cs2 = cs2 * cvel**2
						
					endif
				enddo
			enddo
		enddo
		
		! i = nrho-30
		! j = 1
		! k = nye
		
		! h = 10.0d0**eos_table(i,j,k,ipress) + (10.0d0**eos_table(i,j,k,ienergy) - energy_shift - &
		! 0.511d0 / rmu * ergmev * ye(k) + cvel**2) * 10.0d0**logrho(i)
		
		! kappa_t = 10.0d0**logrho(i) * dpdrho(i,j,k) * (cm3fm3/ergmev)
		! kappa_t = 1.0d0 / kappa_t
		! beta_v = dpdT(i,j,k) * (cm3fm3/ergmev) / kmev
		! cv = 1.0d0 / (10.0d0**logrho(i)/rmu*cm3fm3) * dsdlnT(i,j,k) / LOG(10.0d0) 
		! cp = cv + 10.0d0**logtemp(j)*kmev/(10.0d0**logrho(i)/rmu*cm3fm3) * beta_v * kappa_t * beta_v
		! gamma_ad = cp/cv
						
		! cs2(i,j,k) = gamma_ad / (h * (cm3fm3/ergmev) * kappa_t )
		! gamma(i,j,k) = gamma_ad / 10.0d0**eos_table(i,j,k,ipress)*ergmev/cm3fm3 / kappa_t		
						
		! write(*,*) kappa_t, beta_v, cv, cp, gamma_ad, &
		! 10.0d0**logtemp(j)*kmev/(10.0d0**logrho(i)/rmu*cm3fm3) * beta_v * kappa_t * beta_v
		! write(*,*) cs2(i,j,k)
		! write(*,*) 10.0d0**logtemp(j), dpdT(i,j,k), 10.0d0**eos_table(i,j,k,ipress), (cm3fm3/ergmev) / kmev
		! stop
		
		WRITE(*,*) MAXVAL(cs2), MINVAL(cs2)
		WRITE(*,*) MAXVAL(gamma), MINVAL(gamma)
		WRITE(*,*) MAXVAL(dpdrhoT), MAXVAL(dpdt), MAXVAL(dedt) 
		WRITE(*,*) MINVAL(dpdrhoT), MINVAL(dpdt), MINVAL(dedt) 
		
		DEALLOCATE(dedT)
		DEALLOCATE(dsdlnT)
		DEALLOCATE(dsdlnrho)
		DEALLOCATE(dpdrhoT)
		DEALLOCATE(dedrhoT)
		DEALLOCATE(dpdT)
		DEALLOCATE(dpdrho)
		DEALLOCATE(dpde)
		
	end subroutine derivatives_production

END MODULE wlGammaSoundSpeed