# The directories containing the source files, separated by ':'
# When GNU Make cannot find a prerequisite it will search the directories in the VPATH list (from first to last) and stop at the first directory in which it finds the missing prerequisite.
#-------------------------
vpath %.f90 $(WEAKLIB_DIR)/External/Utilities/UnitTests/GreyVariables:$(WEAKLIB_DIR)/External/Utilities/Opacities/ElectronCaptureEmission:$(WEAKLIB_DIR)/Distributions/OpacitySource:$(WEAKLIB_DIR)/Distributions/Library:$(WEAKLIB_DIR)/Distributions/EOSSource:$(WEAKLIB_DIR)/External/SharedModules
#========================= 
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
# The source files: regardless of where they reside in the source tree,
# VPATH will locate them...
#---------------------------------------------------

BASE_MODS = wlKindModule.o glaquad.o GaussianQuadrature.o fexp.o wlEnergyGridModule.o wlGridModule.o wlThermoStateModule.o wlDependentVariablesModule.o wlExtNumericalModule.o wlExtPhysicalConstantsModule.o wlOpacityFieldsModule.o

EOS_MODS = wlEquationOfStateTableModule.o

PHY_MODS = B85.o

Grey_MODS = GreyVariables.o

BaseIO_MODS = wlIOModuleHDF.o

EOSIO_MODS = wlEOSIOModuleHDF.o

Opacity_MODS = wlOpacityTableModule.o

OpacityIO_MODS = wlOpacityTableIOModuleHDF.o

Opacity_OBJ =  $(BASE_MODS) $(EOS_MODS) $(PHY_MODS) $(BaseIO_MODS) $(EOSIO_MODS) $(Grey_MODS) $(Opacity_MODS) $(OpacityIO_MODS) wlGreyVariablesTEST.o
#===================================================
# Build a Dependency list and an Object list, by replacing the .cpp
# extension to .d for dependency files, and .o for object files.
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

#-------------------
#test: $(Opacity_MODS)
test: $(Opacity_OBJ)
	$(LDRL) $(LDRFLAGS) -o wlGREYTEST $(LIBS) $(Opacity_OBJ)

#===================

#---------------------------------
glaquad.o: wlKindModule.o
GaussianQuadrature.o: wlKindModule.o
wlEnergyGridModule.o: wlKindModule.o
wlGridModule.o: wlKindModule.o
wlThermoStateModule.o: wlKindModule.o
wlDependentVariablesModule.o: wlKindModule.o
wlExtNumericalModule.o: wlKindModule.o
wlExtPhysicalConstantsModule.o: wlKindModule.o wlExtNumericalModule.o
B85.o: wlKindModule.o wlExtPhysicalConstantsModule.o fexp.o
wlOpacityFieldsModule.o: wlKindModule.o

$(EOS_MODS): $(BASE_MODS)
$(Grey_MODS):$(BASE_MODS) $(PHY_MODS)
$(BaseIO_MODS): $(BASE_MODS)
$(EOSIO_MODS): $(EOS_MODS) $(BaseIO_MODS) $(BASE_MODS)
$(Opacity_MODS): $(EOSIO_MODS) $(BaseIO_MODS) $(EOS_MODS) $(BASE_MODS) 
$(OpacityIO_MODS): $(Opacity_MODS) $(EOSIO_MODS) $(BASE_MODS)
wlCreateOpacityTable.o: $(OpacityIO_MODS) $(Opacity_MODS) $(BaseIO_MODS) $(Grey_MODS) $(BASE_MODS)

#==================================
scrub:
	rm -f core *.o *.f90_o *.f_o *.mod *.lst *.cub *.ptx *.i *.T *.diag *__genmod.* *.oo 

clean:
	rm -f *.o *.mod *.lst

