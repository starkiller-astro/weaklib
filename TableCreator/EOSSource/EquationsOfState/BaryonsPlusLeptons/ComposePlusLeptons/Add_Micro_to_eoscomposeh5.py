import numpy as np
import sys
import h5py

file_to_modify     = 'Executables/SFHo_no_ele/eoscompose.h5'
file_to_copy_micro = 'Executables/SFHo/eoscompose.h5'

# you find this in the specific eos.pdf file from the CompOSE database
# each index refers to a different microscopic property
dict_micro = { 'Neutron Effective Mass' : 10041, \
               'Proton Effective Mass' : 11041, \
               'Neutron Self Energy' : 10051, \
               'Proton Self Energy' : 11051 } # this is in eos.pdf

# start reading the table
with h5py.File(file_to_copy_micro, 'r') as f:

  index_micro = f['Micro_qty/index_micro'][:]
  micro       = f['Micro_qty/micro'][:]
  pointsmicro = f['Micro_qty'].attrs['pointsmicro']

with h5py.File(file_to_modify, 'r+') as g:

  del g['Micro_qty/index_micro']
  del g['Micro_qty/micro']

  g.create_dataset('Micro_qty/index_micro', data=index_micro)
  g.create_dataset('Micro_qty/micro', data=micro)
  data = g['Micro_qty'].attrs['pointsmicro']
  data = pointsmicro


