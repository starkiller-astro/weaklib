# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 15:10:55 2024

@author: Lucab
"""

import numpy as np
import sys

def read_compo(file_compo, file_yq, dict_pairs, charge_tol=1e-7, mfrac_tol=1e-7):

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
                if abs(Xh[iL] + Xa[iL] + Xp[iL] + Xn[iL] - 1.0) > mfrac_tol:
                    print('mass fraction not conserved 1', iL, Xh[iL] + Xp[iL] + Xn[iL])
                    sys.exit()
                    
                # check charge conservation            
                if abs(Xp[iL] + Xa[iL]/2. - yq[iyq[iL]-1]) > charge_tol:
                    print('charge not conserved 2', iL, Xp[iL] + Xa[iL]/2., yq[iyq[iL]-1])
                    sys.exit()
                
                continue            
            else:
                index = data[5 + 2*Npairs+1]
                if index != dict_pairs['index_avg_nucleus']:
                    print('I do not know what particle this is')
                    sys.exit()
                    
                Abar[iL] = data[5 + 2*Npairs+2]
                Zbar[iL] = data[5 + 2*Npairs+3]
                Xh[iL] = data[5 + 2*Npairs+4] * Abar[iL]
                
                # check mass fraction sums to 1
                if abs(Xh[iL] + Xa[iL] + Xp[iL] + Xn[iL] - 1.0) > mfrac_tol:
                    print('mass fraction not conserved 3', iL, Xh[iL] + Xp[iL] + Xn[iL])
                    sys.exit()
                    
                # check charge conservation
                if abs(Xp[iL] + Xa[iL]/2. + Zbar[iL]/Abar[iL]*Xh[iL] - yq[iyq[iL]-1]) > charge_tol:
                    print('charge not conserved 3', iL, Xp[iL] + Xa[iL]/2. + Zbar[iL]/Abar[iL]*Xh[iL], yq[iyq[iL]-1])
                    sys.exit()        

    return iT, irho, iyq, Xp, Xn, Xa, Xh, Abar, Zbar

def read_micro(file_micro, dict_micro):
    with open(file_micro, 'r') as f:
        lines = f.readlines()
        
        N = len(lines)
        
        irho, iT, iyq = np.zeros(N,dtype=int), np.zeros(N,dtype=int), np.zeros(N,dtype=int)
        p_self_ene, n_self_ene = np.zeros(N), np.zeros(N)
        for iL, line in enumerate(lines):
            
            elements = [l.strip() for l in line.split()]
            data = [int(ele) if ele.isdigit() else float(ele) for ele in elements]        
            iT[iL], irho[iL], iyq[iL] = data[0], data[1], data[2]

            nMicro = data[3]

            # now get indices and values
            indices = np.array(data[4:][::2])
            values  = np.array(data[4:][1::2])

            in_self_ene = np.where( indices == dict_micro['Neutron Self Energy'] )[0][0]
            n_self_ene[iL] = values[in_self_ene]

            ip_self_ene = np.where( indices == dict_micro['Proton Self Energy'] )[0][0]
            p_self_ene[iL] = values[ip_self_ene]

    return n_self_ene, p_self_ene

# set tolerances for charge and mass conservation
charge_tol = 1e-7
mfrac_tol = 1e-7

file_compo = 'Executables/SFHo_no_ele/eos.compo'
file_micro = 'Executables/SFHo/eos.micro'
file_yq = 'Executables/SFHo_no_ele/eos.yq'
file_WL_format = 'Executables/SFHo_no_ele/eos.compomicro.wl'

# you find this in the specific eos.pdf file from the CompOSE database
# 10 and 11 are always neutron and proton. Each element of the dictionary is:
# 'index' : (charge, baryon number). Technically you can infer from the index 
# the charge and baryon numbers, but then if there are mesons it gets tricky I think.
dict_pairs = { '10'  : (0,1), \
               '11'  : (1,1), \
               '0'   : (1,0), \
               '4002': (2,4), \
               '3002': (2,3), \
               '3001': (1,3), \
               '2001': (1,2), \
               'index_avg_nucleus' : 999} # this is in eos.pdf

# now the microscopic part
dict_micro = { 'Neutron Effective Mass' : 10041, \
               'Proton Effective Mass' : 11041, \
               'Neutron Self Energy' : 10051, \
               'Proton Self Energy' : 11051 } # this is in eos.pdf

n_self_ene, p_self_ene = read_micro(file_micro, dict_micro)

iT, irho, iyq, Xp, Xn, Xa, Xh, Abar, Zbar = read_compo(file_compo, file_yq, dict_pairs, \
                                                      charge_tol=charge_tol, mfrac_tol=mfrac_tol)


data_to_write = np.column_stack((iT, irho, iyq, Xp, Xn, Xa, Xh, Abar, Zbar, n_self_ene, p_self_ene))
np.savetxt(file_WL_format, data_to_write, fmt=['%3.0i ','%3.0i ','%3.0i ', \
    '%.7E ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E '])