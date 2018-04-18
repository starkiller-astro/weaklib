
# Build configurations using variables defined in System Makefile
CMODE = OPTIMIZE
#CMODE = DEBUG
#CMODE = PROFILE

FFLAGS   = $(R8)
CCFLAGS  = 
CXXFLAGS =
LDRFLAGS = $(R8)

# Choose matrix solver package 
MATRIX_SOLVER = dense
#MATRIX_SOLVER = PARDISO
#MATRIX_SOLVER = gpu

# Choose executable name
EXE = gettermeta

#SELECT MPI MODE, ON=MPI, OFF= noMPI 
MPI_MODE = ON
#MPI_MODE ?= OFF

#SELECT OPENMP MODE, Default=OFF
#OPENMP_MODE = ON
OPENMP_MODE ?= OFF

#SELECT GPU MODE, Default=OFF
#GPU_MODE = ON
GPU_MODE ?= OFF

# System configurations
# Architecture and installation specific MACROS
include ./Makefile.internal

ifeq ($(CMODE),OPTIMIZE)
  FFLAGS   += $(OPT) $(EXTRA_FFLAGS)
  CCFLAGS  += $(OPT) $(EXTRA_CCFLAGS)
  CXXFLAGS += $(OPT) $(EXTRA_CXXFLAGS)
  LDRFLAGS += $(OPT) $(EXTRA_LDRFLAGS)
else ifeq ($(CMODE),DEBUG)
  FFLAGS   += $(DEBUG) $(EXTRA_FFLAGS)
  CCFLAGS  += $(DEBUG) $(EXTRA_CCFLAGS)
  CXXFLAGS += $(DEBUG) $(EXTRA_CXXFLAGS)
  LDRFLAGS += $(DEBUG) $(EXTRA_LDRFLAGS)
else ifeq ($(CMODE),PROFILE)
  FFLAGS   += $(PROFILE) $(EXTRA_FFLAGS)
  CCFLAGS  += $(PROFILE) $(EXTRA_CCFLAGS)
  CXXFLAGS += $(PROFILE) $(EXTRA_CXXFLAGS)
  LDRFLAGS += $(PROFILE) $(EXTRA_LDRFLAGS)
endif

MODULE_SRC = KIND_MODULE.o \
             NUMERICAL_MODULE.f90 \
             PARALLEL_MODULE.f90 \
             PHYSCNST_MODULE.f90 \
             PRB_CNTL_MODULE.f90 \
             RADIAL_RAY_MODULE.f90 \
             ARRAY_MODULE.f90 \
             CYCLE_MODULE.f90 \
             EDIT_MODULE.f90 \
             NU_DIST_MODULE.f90 \
             NU_ENERGY_GRID_MODULE.f90 \
             EL_EOS_MODULE.f90 \
             EOS_BCK_MODULE.f90 \
             EOS_LS_MODULE.f90 \
             EOS_M4C_MODULE.f90 \
             EOS_SNC_X_MODULE.f90 \
             EVH1_SWEEP_MODULE.f90 \
             MAXWEL_MODULE.f90 \
             NUCBRN_MODULE.f90 \
             TEST_MODULE.f90 \
             EC_TABLE_MODULE.f90 \
             POLYLOG_MODULE.f90 \
             deflash_x.f90 \
             e_grid_interpolation.f90 \
             vector_functions.f90 \
             dimension_edit_arrays.f90 \
             dimension_eos_bck_arrays.f90 \
             dimension_eos_ls_arrays.f90 \
             dimension_eos_snc_x_arrays.f90 \
             dimension_nucbrn_arrays.f90 \
             Eos_bck.f \
             eos.f90 \
             eosidx_x.f90 \
             Eos_ls.f \
             eoslsgen_x.f90 \
             eos_m_gamma.f90 \
             eosnuc_x.f90 \
             eos_x_gamma.f90 \
             e_p_eos.f90 \
             eqstlsa_x.f90 \
             eqstt_x.f90 \
             eqstz_x.f90 \
             esrgn_x.f90 \
             esrgn_x_t.f90 \
             esrgnz_x.f90 \
             fexp.f90 \
             flush_xt.f90 \
             glaquad.f90 \
             lectron.f90 \
             lectron_gamma.f90 \
             loadmx.f \
             net_comp.f90 \
             net.f90 \
             numerical_assign.f90 \
             terminate.f90 \
             tgvndeye_x.f90 \
             abem_cal.f90 \
             abemfrnpetanp.f90 \
             abemfrnp.f90 \
             ABEM_MODULE.o \
             abemnc.f90 \
             abemrate.f90 \
             abemrgn.f90 \
             abemset.f90 \
             abemtnc.f90 \
             cc_weak_mag.f90 \
             interp_ec.f90 \
             nu_N_absr_momts.f90 \
             f10.f90 \
             gquad.f90 \
             index_transpose.f90 \
             cc_difcs.f \
             n_difcs.f \
             load_polylog.f90 \
             EC_TABLE_MODULE.o 

EXEC_SRC = eos_get.f90 table_get.f90

ALL_SRC = $(MODULE_SRC) $(EXEC_SRC) $(EXTRA_SRC)
ALL_OBJ := $(ALL_SRC:.f90=.o)
ALL_OBJ := $(ALL_OBJ:.f=.o)  

# Configure Libraries and files to link Matrix solver
LIBS = $(EXTRA_LINK)


ifeq ($(MATRIX_SOLVER),dense)
  LIBS += $(LAPACK_LIBDIR) $(LAPACK_LIBS)
  EXTRA_OBJ += $(LAPACK_OBJ)

else ifeq ($(MATRIX_SOLVER),PARDISO)
  LIBS += $(PARDISO_LIBDIR) $(PARDISO_LIBS)
  EXTRA_OBJ += $(PARDISO_OBJ)

else ifeq ($(MATRIX_SOLVER),gpu)
  LIBS += $(LAPACK_LIBDIR) $(LAPACK_LIBS)
  EXTRA_OBJ += $(LAPACK_OBJ)
  LIBS       += $(GPU_LAPACK_LIBDIR) $(GPU_LAPACK_LIBS)
  EXTRA_OBJ  += $(GPU_LAPACK_OBJ)
  LAPACK_INC += $(GPU_LAPACK_INC)

endif

# Configure flags for OpenMP
FFLAGS   += $(OMP_FLAG)
CCFLAGS  += $(OMP_FLAG)
CXXFLAGS += $(OMP_FLAG)
LDRFLAGS += $(OMP_FLAG)

# Configure flags for OpenACC
ifeq ($(GPU_MODE),ON)
  ifeq ($(MATRIX_SOLVER),dense)
    LIBS       += $(GPU_LAPACK_LIBDIR) $(GPU_LAPACK_LIBS)
    EXTRA_OBJ  += $(GPU_LAPACK_OBJ)
    LAPACK_INC += $(GPU_LAPACK_INC)
  endif
  ifeq ($(PE_ENV),PGI)
    ifeq ($(ACCMODE),OPENACC)
      CCFLAGS  += -acc -ta=nvidia
      CXXFLAGS += -acc -ta=nvidia
      FFLAGS   += -acc -ta=nvidia
      LDRFLAGS += -acc -ta=nvidia
    else
      CCFLAGS  += -ta=nvidia
      CXXFLAGS += -ta=nvidia
      FFLAGS   += -ta=nvidia
      LDRFLAGS += -ta=nvidia
    endif
  endif
else
  ifeq ($(PE_ENV),CRAY)
    CCFLAGS  += -h noacc
    CXXFLAGS += -h noacc
    FFLAGS   += -h noacc
    LDRFLAGS += -h noacc
  endif
endif

FCL    = $(FC)
CCL    = $(CC)
CXXL   = $(CXX)
LDRL   = $(LDR)

$(EXE): $(ALL_OBJ)
	$(LDRL) $(LDRFLAGS) -o $(EXE) $(ALL_OBJ) $(LIBS)
#
# Rules for compiling individual files.
#
$($(LAPACK_VER)_OBJ): %.o: %.f90
	$(FC) $(FFLAGS) $($(LAPACK_VER)_INC) \
	    -c $< -o $(notdir $@)
$($(GPU_LAPACK_VER)_OBJ): %.o: %.f90
	$(FC) $(FFLAGS) $($(GPU_LAPACK_VER)_INC) \
	    -c $< -o $(notdir $@)
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(CUDATK_INC)
%.o: %.f
	$(FCL)  $(FFLAGS) -c $< -o $@
%.o: %.f90
	$(FCL) $(FFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: wlThermoStateModule.o wlThermoStateTest.o wlKindModule.o 
	$(LDRL) $(LDRFLAGS) -o wltest $(LIBS) wlKindModule.o wlThermoStateModule.o wlThermoStateTest.o

clean:
	rm -f core *.o *.f90_o *.f_o *.mod *.lst *.cub *.ptx *.i *.T *.diag *__genmod.* *.oo 

