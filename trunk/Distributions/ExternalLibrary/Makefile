MACHINE = $(WEAKLIB_MACHINE)

OPT_LEVEL = OPTIMIZE
FLAGS     = $(FLAGS_$(OPT_LEVEL))

include ../Build/Makefile_Build

ExternalLibrary: $(weaklib)

clean:
	rm -f *.o *.mod *.a
