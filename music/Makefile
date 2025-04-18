##############################################################################
### compile time configuration options
FFTW3		= yes
MULTITHREADFFTW	= yes
SINGLEPRECISION	= no
HAVEHDF5        = yes
HAVEBOXLIB	= no
BOXLIB_HOME     = ${HOME}/nyx_tot_sterben/BoxLib
FFTW_INCL= -I/user1/poojarani/fftw-2.1.5/include
FFTW_LIBS= -L/user1/poojarani/fftw-2.1.5/lib
GSL_INCL= -I/user1/poojarani/gsl-2.5/include
GSL_LIBS= -L/user1/poojarani/gsl-2.5/lib
HDF5_INCL= -I/user1/poojarani/hdf5-1.13.1/include
HDF5_LIBS= -L/user1/poojarani/hdf5-1.13.1/lib
##############################################################################
### compiler and path settings
CC      = g++ -std=c++11
OPT     = -Wall -Wno-unknown-pragmas -O3 -g -mtune=native
CFLAGS  =  
LFLAGS  = -lgsl -lgslcblas 
CPATHS  = -I./src -I$(HOME)/local/include -I/opt/local/include -I/usr/local/include $(FFTW_INCL) $(GSL_INCL) $(HDF5_INCL)
LPATHS  = -L$(HOME)/local/lib -L/opt/local/lib -L/usr/local/lib $(FFTW_LIBS) $(GSL_LIBS) $(HDF5_LIBS)

##############################################################################
# if you have FFTW 2.1.5 or 3.x with multi-thread support, you can enable the 
# option MULTITHREADFFTW
ifeq ($(strip $(MULTITHREADFFTW)), yes)
  ifeq ($(CC), mpiicpc)
    CFLAGS += -openmp
    LFLAGS += -openmp
  else
    CFLAGS += -fopenmp
    LFLAGS += -fopenmp
  endif
  ifeq ($(strip $(FFTW3)),yes)
	ifeq ($(strip $(SINGLEPRECISION)), yes)
		LFLAGS  +=  -lfftw3f_threads
	else
		LFLAGS  +=  -lfftw3_threads
	endif
  else
    ifeq ($(strip $(SINGLEPRECISION)), yes)
      LFLAGS  += -lsrfftw_threads -lsfftw_threads
    else
      LFLAGS  += -ldrfftw_threads -ldfftw_threads
    endif
  endif
else
  CFLAGS  += -DSINGLETHREAD_FFTW
endif

ifeq ($(strip $(FFTW3)),yes)
  CFLAGS += -DFFTW3
endif

##############################################################################
# this section makes sure that the correct FFTW libraries are linked
ifeq ($(strip $(SINGLEPRECISION)), yes)
  CFLAGS  += -DSINGLE_PRECISION
  ifeq ($(FFTW3),yes)
    LFLAGS += -lfftw3f
  else
    LFLAGS  += -lsrfftw -lsfftw
  endif
else
  ifeq ($(strip $(FFTW3)),yes)
    LFLAGS += -lfftw3
  else
    LFLAGS  += -ldrfftw -ldfftw
  endif
endif

##############################################################################
#if you have HDF5 installed, you can also enable the following options
ifeq ($(strip $(HAVEHDF5)), yes)
  OPT += -DH5_USE_16_API -DHAVE_HDF5
  LFLAGS += -lhdf5
endif

##############################################################################
CFLAGS += $(OPT)
TARGET  = MUSIC
OBJS    = output.o transfer_function.o Numerics.o defaults.o constraints.o random.o\
		convolution_kernel.o region_generator.o densities.o cosmology.o poisson.o\
		densities.o cosmology.o poisson.o log.o main.o \
		$(patsubst src/plugins/%.cc,src/plugins/%.o,$(wildcard src/plugins/*.cc))

##############################################################################
# stuff for BoxLib
BLOBJS = ""
ifeq ($(strip $(HAVEBOXLIB)), yes)
  IN_MUSIC = YES
  TOP = ${PWD}/src/plugins/nyx_plugin
  CCbla := $(CC)
  include src/plugins/nyx_plugin/Make.ic
  CC  := $(CCbla)
  CPATHS += $(INCLUDE_LOCATIONS)
  LPATHS += -L$(objEXETempDir)
  BLOBJS = $(foreach obj,$(objForExecs),src/plugins/boxlib_stuff/$(obj))
#
endif

##############################################################################
all: $(OBJS) $(TARGET) Makefile
#	cd plugins/boxlib_stuff; make

bla:
	echo $(BLOBJS)

blabla:
	echo $(OBJS)

ifeq ($(strip $(HAVEBOXLIB)), yes)
$(TARGET): $(OBJS) src/plugins/nyx_plugin/*.cpp
	cd src/plugins/nyx_plugin; make BOXLIB_HOME=$(BOXLIB_HOME) FFTW3=$(FFTW3) SINGLE=$(SINGLEPRECISION)
	$(CC) $(LPATHS) -o $@ $^ $(LFLAGS) $(BLOBJS) -lifcore
else
$(TARGET): $(OBJS)
	$(CC) $(LPATHS) -o $@ $^ $(LFLAGS)
endif

%.o: src/%.cc src/*.hh Makefile 
	$(CC) $(CFLAGS) $(CPATHS) -c $< -o $@

src/plugins/%.o: src/plugins/%.cc src/*.hh Makefile 
	$(CC) $(CFLAGS) $(CPATHS) -c $< -o $@

clean:
	rm -rf $(OBJS)
ifeq ($(strip $(HAVEBOXLIB)), yes)
	oldpath=`pwd`
	cd src/plugins/nyx_plugin; make realclean BOXLIB_HOME=$(BOXLIB_HOME)
endif
	cd $(oldpath)
	
