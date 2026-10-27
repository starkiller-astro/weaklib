MACHINE = $(WEAKLIB_MACHINE)

OPT_LEVEL = OPTIMIZE
FLAGS     = $(FLAGS_$(OPT_LEVEL))

USE_GPU    ?= FALSE
USE_CUDA   ?= FALSE
USE_HIP    ?= FALSE
USE_OMP_OL ?= FALSE
USE_OACC   ?= FALSE
USE_OMP    ?= TRUE
USE_CUBA   ?= FALSE

# COMPOSE: the program reads BaryonsPlusPhotonsPlusLeptonsEOS.h5 and needs the
# muon chemical potential (in 3D mode ApplyEOS sets Mumu = 0)
EOSMODE = COMPOSE

include ../Makefile_CreateTable
