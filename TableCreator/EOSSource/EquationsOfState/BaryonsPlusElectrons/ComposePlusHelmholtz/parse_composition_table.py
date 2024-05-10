# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 15:10:55 2024

@author: Lucab
"""

import numpy as np
import sys

# set tolerances for charge and mass conservation
charge_tol = 1e-7
mass_frac_tol = 1e-7

Compose_dir = 'SFHo'
Compose_dir = 'SFHo_no_ele'

file_compo = Compose_dir + '/eos.compo'
file_yq = Compose_dir + '/eos.yq'
file_for_WL = Compose_dir + '/eos.compo.wl'

# you find this in the specific eos.pdf file from the CompOSE database
# 10 and 11 are always neutron and proton. Each element of the dictionary is:
# 'index' : (charge, baryon number). Technically you can infer from the index 
# the charge and baryon numbers, but then if there are mesons it gets tricky I think.
dict_pairs = { '10': (0,1), \
               '11': (1,1), \
               '0': (1,0), \
               '4002': (2,4), \
               '3002': (2,3), \
               '3001': (1,3), \
               '2001': (1,2) }

index_avg_nucleus = 999 # this is in eos.pdf
    
yq = np.loadtxt(file_yq, skiprows=2)

with open(file_compo, 'r') as f:
    lines = f.readlines()
    
    N = len(lines)
    
    irho, iT, iyq = np.zeros(N,dtype=int), np.zeros(N,dtype=int), np.zeros(N,dtype=int)
    Xp, Xn, Xa, Xh, Abar, Zbar = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
    for iL, line in enumerate(lines):
        
        elements = [l.strip() for l in line.split()]
        data = [int(ele) if ele.isdigit() else float(ele) for ele in elements]
        
        # now that you have the data, identify how many pairs and quads you have
        iT[iL], irho[iL], iyq[iL] = data[0], data[1], data[2]
        iphase, Npairs = data[3], data[4]
        
        if len(data) > 5 + 2*Npairs:
            Nquad = data[5 + 2*Npairs]
        else:
            Nquad = 0
            if len(data) != 5 + 2*Npairs:
                print('incompatible data length 1', len(data), 5 + 2*Npairs )
                sys.exit()
                
        if len(data) != 5 + 2*Npairs + 4*Nquad+1:
            print('incompatible data length 2', len(data), 5 + 2*Npairs + 4*Nquad+1 )
            sys.exit()
            
        for ipair in range(Npairs):
            dict_index = data[5 + 2*ipair]
            abundance_fraction = data[5 + 2*ipair+1]
            baryon_number = dict_pairs[str(dict_index)][1]
            mass_fraction = abundance_fraction * baryon_number
            charge = dict_pairs[str(dict_index)][0]
            
            if dict_index == 10:
                Xn[iL] = mass_fraction
            elif dict_index == 11:
                Xp[iL] = mass_fraction
            elif dict_index == 4002:
                Xa[iL] = mass_fraction
            elif dict_index == 0:
                electron_fraction = abundance_fraction
            else:
                # lumping light nuclei in to protons and neutrons
                Xp[iL] += mass_fraction * charge / baryon_number
                Xn[iL] += mass_fraction * (baryon_number - charge) / baryon_number
        
        if Nquad > 1:
            print("I can't handle this, sorry")
            # sys.exit()
        
        if Nquad == 0:
            # check mass fraction sums to 1
            if abs(Xh[iL] + Xa[iL] + Xp[iL] + Xn[iL] - 1.0) > mass_frac_tol:
                print('mass fraction not conserved 1', iL, Xh[iL] + Xp[iL] + Xn[iL])
                sys.exit()
                
            # check charge conservation            
            if abs(Xp[iL] + Xa[iL]/2. - yq[iyq[iL]-1]) > charge_tol:
                print('charge not conserved 2', iL, Xp[iL] + Xa[iL]/2., yq[iyq[iL]-1])
                sys.exit()
            
            continue            
        else:
            index = data[5 + 2*Npairs+1]
            if index != index_avg_nucleus:
                print('I do not know what particle this is')
                sys.exit()
                
            Abar[iL] = data[5 + 2*Npairs+2]
            Zbar[iL] = data[5 + 2*Npairs+3]
            Xh[iL] = data[5 + 2*Npairs+4] * Abar[iL]
            
            # check mass fraction sums to 1
            if abs(Xh[iL] + Xa[iL] + Xp[iL] + Xn[iL] - 1.0) > mass_frac_tol:
                print('mass fraction not conserved 3', iL, Xh[iL] + Xp[iL] + Xn[iL])
                sys.exit()
                
            # check charge conservation
            if abs(Xp[iL] + Xa[iL]/2. + Zbar[iL]/Abar[iL]*Xh[iL] - yq[iyq[iL]-1]) > charge_tol:
                print('charge not conserved 3', iL, Xp[iL] + Xa[iL]/2. + Zbar[iL]/Abar[iL]*Xh[iL], yq[iyq[iL]-1])
                sys.exit()        

np.savetxt(file_for_WL, np.column_stack((iT, irho, iyq, Xp, Xn, Xa, Xh, Abar, Zbar)), \
           fmt=['%3.0i ','%3.0i ','%3.0i ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E '])