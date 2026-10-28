MACHINE = $(WEAKLIB_MACHINE)

OPT_LEVEL = OPTIMIZE
FLAGS     = $(FLAGS_$(OPT_LEVEL))

USE_GPU    ?= FALSE
USE_CUDA   ?= FALSE
USE_HIP    ?= FALSE
USE_OMP_OL ?= FALSE
USE_OACC   ?= FALSE
USE_OMP    ?= TRUE

# COMPOSE: wlTestElasticOpacity (reads BaryonsPlusPhotonsPlusLeptonsEOS.h5).
# wlTestGeneralKernelVsTable(_Interface) need EOSMODE=3D (SFHo EOS):
#   make clobber; make wlTestGeneralKernelVsTable EOSMODE=3D
EOSMODE ?= COMPOSE

include ../Makefile_GeneralKernelsTest
