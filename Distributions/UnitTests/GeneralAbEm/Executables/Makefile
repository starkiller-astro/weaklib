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

# 3D: old weaklib SFHo table (needed by wlTestAbEm2DVsTable)
# COMPOSE: BaryonsPlusPhotonsPlusLeptonsEOS.h5 (compare targets)
# Run 'make clobber' when switching EOSMODE or USE_CUBA.
EOSMODE ?= 3D

include ../Makefile_GeneralAbEm
