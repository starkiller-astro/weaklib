# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 15:10:55 2024

@author: Lucab
"""

import numpy as np
import sys

DummyValue = 123456

def HandleLightClusters(XDict, ZDict, ADict, dict_compo_names, How = 'LumpIntoNucleons'):
    
    if How == 'LumpIntoNucleons':
        return LumpIntoNucleons(XDict, ZDict, ADict, dict_compo_names)
    elif How == 'LumpIntoHeavyNucleus':
        return LumpIntoHeavyNucleus(XDict, ZDict, ADict, dict_compo_names)
    elif How == 'LumpIntoAlphasAndNucleons':
        return LumpIntoAlphasAndNucleons(XDict, ZDict, ADict, dict_compo_names)
    else:
        raise RuntimeError('Specify valid strtegy to handle light nuclei')
  
def LumpIntoNucleons(XDict, ZDict, ADict, dict_compo_names):
    
    Xp = XDict[dict_compo_names['Proton']]
    Xn = XDict[dict_compo_names['Neutron']]
    Xa = XDict[dict_compo_names['Alpha']]

    indices = list(XDict.keys())
    names = list(dict_compo_names.keys())
    ClustersIndices = [dict_compo_names[name] for name in names \
                       if not name in ['Proton', 'Neutron', 'Alpha']]

    for key in ClustersIndices:
        mass_fraction, charge, baryon_number = XDict[key], ZDict[key], ADict[key]

        Xp += mass_fraction * charge / baryon_number
        Xn += mass_fraction * (baryon_number - charge) / baryon_number

    return XDict['Heavy'], ZDict['Heavy'], ADict['Heavy'], Xp, Xn, Xa

def LumpIntoAlphasAndNucleons(XDict, ZDict, ADict, dict_compo_names):
    # Still not implemented, maybe put one neutron from H3 and/or deuteron into He3 and
    # the rest lump into nucleons. But what if there is a ton of He3 and no deuteron/H3?
    return

def LumpIntoHeavyNucleus(XDict, ZDict, ADict, dict_compo_names):
    
    indices = list(XDict.keys())
    Xp = XDict[dict_compo_names['Proton']]
    Xn = XDict[dict_compo_names['Neutron']]
    Xa = XDict[dict_compo_names['Alpha']]

    indices = list(XDict.keys())
    names = list(dict_compo_names.keys())
    ClustersIndices = [dict_compo_names[name] for name in names \
                       if not name in ['Proton', 'Neutron', 'Alpha']]

    if not 'Heavy' in indices:
        raise RuntimeError('No Heavies it seems in this line')
    
    Xh   = XDict['Heavy']
    Zbar = ZDict['Heavy']
    Abar = ADict['Heavy']

    Zbar_new, Abar_new, Xh_new = 0, 0, 0.0
    for key in ClustersIndices:
        mass_fraction, charge, baryon_number = XDict[key], ZDict[key], ADict[key]

        if mass_fraction == 0:
            continue
        
        # If you get past this then you should have Heavies AND Light Clusters
        if Zbar == DummyValue:
            # if you don't have heavy nuclei but only light clusters it might be
            # that the light cluster abundances are super small. Doesn't really
            # matter what you do with them
            XClust = [XDict[ind] for ind in ClustersIndices]
            if sum(XClust) > 1e-10:
                print('No Heavy Nuclei here, and Light clusters have X =',XDict[ClustersIndices])
            return LumpIntoNucleons(XDict, ZDict, ADict, dict_compo_names)

        Xh_new   = Xh + mass_fraction
        Abar_new = (Xh*Abar      + mass_fraction*baryon_number*mass_fraction) / Xh_new
        Zbar_new = (Xh*Zbar/Abar + mass_fraction*charge/baryon_number)        / Xh_new * Abar_new

        Zbar, Abar, Xh = Zbar_new, Abar_new, Xh_new

    return Xh, Zbar, Abar, Xp, Xn, Xa
    
def read_compo(file_compo, file_yq, dict_pairs, dict_compo_names, \
               HowToHandleLightClusters='LumpIntoNucleons', charge_tol=1e-7, mfrac_tol=1e-7):
    
    yq = np.loadtxt(file_yq, skiprows=2)
    with open(file_compo, 'r') as f:
        lines = f.readlines()

    N = len(lines)
    irho = np.zeros(N, dtype=int)
    iT = np.zeros(N, dtype=int)
    iyq = np.zeros(N, dtype=int)
    Xp = np.zeros(N)
    Xn = np.zeros(N)
    Xa = np.zeros(N)
    Xh = np.zeros(N)
    Abar = np.zeros(N)
    Zbar = np.zeros(N)

    for iL, line in enumerate(lines):
        elements = line.split()
        data = [int(ele) if ele.isdigit() else float(ele) for ele in elements]

        iT[iL], irho[iL], iyq[iL] = data[0], data[1], data[2]
        iphase, Npairs = data[3], data[4]

        base_idx = 5
        if len(data) > base_idx + 2 * Npairs:
            Nquad = data[base_idx + 2 * Npairs]
            has_quad = True
        else:
            Nquad = 0
            has_quad = False
            if len(data) != base_idx + 2 * Npairs:
                raise RuntimeError(f'incompatible data length 1 at line {iL}: {len(data)} vs {base_idx + 2 * Npairs}')

        expected_len = base_idx + 2 * Npairs + (1 + 4 * Nquad if has_quad else 0)
        if len(data) != expected_len:
            raise RuntimeError(f'incompatible data length 2 at line {iL}: {len(data)} vs {expected_len}')

        # Initialize dictionaries with composition information
        XDict, ZDict, ADict = {}, {}, {}
        for key in list(dict_pairs.keys())[:-1]:
            XDict[key] = 0.0
            charge, baryon_number = dict_pairs[key]
            ZDict[key] = charge
            ADict[key] = baryon_number

        XDict['Heavy'] = 0.0
        ZDict['Heavy'] = DummyValue # Should not matter!
        ADict['Heavy'] = DummyValue # Should not matter!

        for ipair in range(Npairs):
            dict_index = data[base_idx + 2 * ipair]
            abundance_fraction = data[base_idx + 2 * ipair + 1]

            key = str(dict_index)
            if key not in dict_pairs:
                raise RuntimeError(f"Unknown index {key} at line {iL}")
            
            charge, baryon_number = dict_pairs[key]
            mass_fraction = abundance_fraction * baryon_number

            XDict[key] = mass_fraction
            ZDict[key] = charge
            ADict[key] = baryon_number

        if Nquad > 1:
            raise RuntimeError(f"I can't handle Nquad > 1 at line {iL}")

        if Nquad == 0:
            Xh[iL], _, _, Xp[iL], Xn[iL], Xa[iL] = HandleLightClusters(XDict, ZDict, ADict, \
                                      dict_compo_names, How=HowToHandleLightClusters)
            if Xh[iL] > 0:
                raise RuntimeError('This should not happen, investigate!')
            msum = Xh[iL] + Xa[iL] + Xp[iL] + Xn[iL]
            if abs(msum - 1.0) > mfrac_tol:
                raise RuntimeError(f"mass fraction not conserved at line {iL}: {msum}")

            charge_sum = Xp[iL] + Xa[iL]/2.
            if abs(charge_sum - yq[iyq[iL]-1]) > charge_tol:
                raise RuntimeError(f"charge not conserved at line {iL}: {charge_sum} vs {yq[iyq[iL]-1]}")
            continue

        quad_idx = base_idx + 2 * Npairs + 1
        index = data[quad_idx]
        if index != dict_pairs['index_avg_nucleus']:
            raise RuntimeError(f"Unknown quad index {index} at line {iL}")

        Abar[iL] = data[quad_idx + 1]
        Zbar[iL] = data[quad_idx + 2]
        Xh[iL]   = data[quad_idx + 3] * Abar[iL]

        XDict['Heavy'] = Xh[iL]
        ZDict['Heavy'] = Zbar[iL]
        ADict['Heavy'] = Abar[iL]
        
        Xh[iL], Zbar[iL], Abar[iL], Xp[iL], Xn[iL], Xa[iL] = \
          HandleLightClusters(XDict, ZDict, ADict, \
              dict_compo_names, How=HowToHandleLightClusters)

        msum = Xh[iL] + Xa[iL] + Xp[iL] + Xn[iL]
        if abs(msum - 1.0) > mfrac_tol:
            print(data)
            raise RuntimeError(f"mass fraction not conserved (quad) at line {iL}: {msum}")

        charge_sum = Xp[iL] + Xa[iL]/2. + Zbar[iL]/Abar[iL]*Xh[iL]
        if abs(charge_sum - yq[iyq[iL]-1]) > charge_tol:
            print(data)
            raise RuntimeError(f"charge not conserved (quad) at line {iL}: {charge_sum} vs {yq[iyq[iL]-1]}")

    return iT, irho, iyq, Xp, Xn, Xa, Xh, Abar, Zbar

def read_micro(file_micro, dict_micro):
    
    with open(file_micro, 'r') as f:
        lines = f.readlines()

    N = len(lines)
    irho = np.zeros(N, dtype=int)
    iT = np.zeros(N, dtype=int)
    iyq = np.zeros(N, dtype=int)
    p_eff_mass = np.zeros(N)
    n_eff_mass = np.zeros(N)
    p_self_ene = np.zeros(N)
    n_self_ene = np.zeros(N)

    id_n_eff_mass = dict_micro['Neutron Effective Mass']
    id_p_eff_mass = dict_micro['Proton Effective Mass']
    id_n_self_ene = dict_micro['Neutron Self Energy']
    id_p_self_ene = dict_micro['Proton Self Energy']

    for iL, line in enumerate(lines):
        # Keep original parsing logic
        elements = line.split()
        data = [int(ele) if ele.isdigit() else float(ele) for ele in elements]

        iT[iL], irho[iL], iyq[iL] = data[0], data[1], data[2]
        nMicro = data[3]

        indices = data[4:][::2]
        values  = data[4:][1::2]

        # Use dictionary for fast lookup
        id_val = dict(zip(indices, values))

        try:
            n_eff_mass[iL] = id_val[id_n_eff_mass]
            p_eff_mass[iL] = id_val[id_p_eff_mass]
            n_self_ene[iL] = id_val[id_n_self_ene]
            p_self_ene[iL] = id_val[id_p_self_ene]
        except KeyError as e:
            raise RuntimeError(f"Missing micro value at line {iL}: {e}")

    return n_eff_mass, p_eff_mass, n_self_ene, p_self_ene

# set tolerances for charge and mass conservation
charge_tol = 1e-7
mfrac_tol = 1e-7

file_yq = 'Executables/SFHo_no_ele/eos.yq'
file_compo = 'Executables/SFHo_no_ele/eos.compo'
file_micro = 'Executables/SFHo/eos.micro'
# Sometimes file_micro needs to be the full table 
# since effective mass and nuclear potentials are 
# not always included in the table without electrons

file_WL_format = 'Executables/SFHo_no_ele/eos.compomicro.wl'

# you find this in the specific eos.pdf file from the CompOSE database
# 10 and 11 are always neutron and proton. Each element of the dictionary is:
# 'index' : (charge, baryon number). Technically you can infer from the index 
# the charge and baryon numbers, but then if there are mesons it gets tricky I think.
dict_pairs = { '10'  : (0,1), \
               '11'  : (1,1), \
               '4002': (2,4), \
               '3002': (2,3), \
               '3001': (1,3), \
               '2001': (1,2), \
               'index_avg_nucleus' : 999} # this is in eos.pdf

dict_compo_names = {'Neutron'  : '10'  , \
                    'Proton'   : '11'  , \
                    'Alpha'    : '4002', \
                    'He3'      : '3002', \
                    'H3'       : '3001', \
                    'Deuteron' : '2001'  \
                    }

HowToHandleLightClusters = 'LumpIntoHeavyNucleus'
# HowToHandleLightClusters = 'LumpIntoNucleons'

iT, irho, iyq, Xp, Xn, Xa, Xh, Abar, Zbar = read_compo(file_compo, file_yq, dict_pairs, dict_compo_names, \
                            HowToHandleLightClusters=HowToHandleLightClusters, charge_tol=charge_tol, mfrac_tol=mfrac_tol)

# now the microscopic part
dict_micro = { 'Neutron Effective Mass' : 10041, \
               'Proton Effective Mass' : 11041, \
               'Neutron Self Energy' : 10051, \
               'Proton Self Energy' : 11051 } # this is in eos.pdf

n_eff_mass, p_eff_mass, n_self_ene, p_self_ene = \
  read_micro(file_micro, dict_micro)

data_to_write = np.column_stack((iT, irho, iyq, Xp, Xn, Xa, Xh, Abar, Zbar, \
                                 n_eff_mass, p_eff_mass, n_self_ene, p_self_ene))
np.savetxt(file_WL_format, data_to_write, fmt=['%3.0i ','%3.0i ','%3.0i ', \
    '%.7E ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E ','%.7E '])