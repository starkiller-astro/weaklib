vpath %.f90 ../../LS:../../BCK:../../SharedModules:../../../Distributions/Source:..
vpath %.f ../../LS:../../BCK
# Build configurations using variables defined in System Makefile
#CMODE = OPTIMIZE
CMODE = DEBUG
#CMODE = PROFILE

FFLAGS   = $(R8)
CCFLAGS  = 
CXXFLAGS =
LDRFLAGS = $(R8)

#SELECT MPI MODE, ON=MPI, OFF= noMPI 
#MPI_MODE = ON
MPI_MODE ?= OFF

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

# Configure Libraries and files to link Matrix solver
LIBS = $(EXTRA_LINK) $(HDF5_FLIB)


# Configure flags for OpenMP
FFLAGS   += $(OMP_FLAG)
CCFLAGS  += $(OMP_FLAG)
CXXFLAGS += $(OMP_FLAG)
LDRFLAGS += $(OMP_FLAG)

# Configure flags for OpenACC
ifeq ($(GPU_MODE),ON)
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

#  Groups of objects

BASE_MODS = wlKindModule.o wlExtNumericalModule.o wlExtPhysicalConstantsModule.o fexp.o

EOSInit_OBJ = initferm.o loadmx.o wlExtInitializeEOS.o

EOS_MODS  = EL_EOS_MODULE.o EOS_M4C_MODULE.o MAXWEL_MODULE.o EOS_BCK_MODULE.o
EOS_OBJ  = LattimerSwesty.o e_p_eos.o BCK_EoS.o glaquad.o net.o

Interface_OBJ = wlExtEOSWrapperModule.o

Interrogator_OBJ = $(BASE_MODS) $(EOS_MODS) $(EOS_OBJ) $(EOSInit_OBJ) $(Interface_OBJ) wlInterrogateEOS.o

FCL    = $(FC)
FCL77  = $(FC77)
CCL    = $(CC)
CXXL   = $(CXX)
LDRL   = $(LDR)

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

interrogate: $(Interrogator_OBJ)
	$(LDRL) $(LDRFLAGS) -o interrogator $(LIBS) $(Interrogator_OBJ)

wlExtNumericalModule.o: wlKindModule.o
wlExtPhysicalConstantsModule.o: wlKindModule.o wlExtNumericalModule.o
$(EOSInit_OBJ): $(EOS_MODS) $(BASE_MODS)
$(EOS_OBJ): $(EOS_MODS) $(BASE_MODS)
$(Interface_OBJ): $(EOS_OBJ) $(EOS_MODS) $(BASE_MODS)
wlInterrogateEOS.o: $(Interface_OBJ) $(EOS_OBJ) $(EOS_MODS) $(BASE_MODS)

scrub:
	rm -f core *.o *.f90_o *.f_o *.mod *.lst *.cub *.ptx *.i *.T *.diag *__genmod.* *.oo 

clean:
	rm -f *.o *.mod *.lst

