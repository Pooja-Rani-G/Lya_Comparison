#############################################################################
# Fortran compiler options and directives

#SYSTYPE="CCF"
#SYSTYPE="lyocral"
#SYSTYPE="lyocralP"
#SYSTYPE="mac"
SYSTYPE="macP"
#SYSTYPE="ccage"
#SYSTYPE="horizon"

OPTIONS = -cpp -fopenmp
OPTIONS += -DDISPLAY_PROGRESS_PERCENT


ifeq ($(SYSTYPE),"lyocral")
F90         = ifort
#FFLAGS      = -g -traceback -fpp -check all -debug -warn all -ftrapuv 
##FFLAGS      = -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS      = -O3 -fpp -ftz -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
##FFLAGS      = -p -g
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
endif

ifeq ($(SYSTYPE),"lyocralP")
F90         = mpif90 -ffree-line-length-none -ffree-form
FFLAGS      = #-g -fcheck=all -ftrapv -fbacktrace -Wall
FFLAGS      = -O2
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
endif


ifeq ($(SYSTYPE),"CCF")
F90         = ifort
##FFLAGS      = -g -fpp -check all -debug -warn all -ftrapuv -traceback
FFLAGS      = -O3 -fpp -ftz  -xavx -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS += $(OPTIONS)
LDFLAGS     =  -lmpi 
endif

ifeq ($(SYSTYPE),"mac")
#F90         = gfortran-mp-4.7 -ffree-line-length-none -ffree-form  ###-fopenmp
F90         = gfortran -ffree-line-length-none -ffree-form  ###-fopenmp
FFLAGS      = -g -fcheck=all -ftrapv -fbacktrace -Wall -fbounds-check
#FFLAGS      = -O2
FFLAGS += $(OPTIONS)
LDFLAGS     = #-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/
endif


ifeq ($(SYSTYPE),"macP")
F90         = mpif90#-openmpi-gcc8
FFLAGS      = -g3 -fcheck=all -ffree-line-length-none -ffree-form -ftrapv -fbounds-check -fbacktrace -Wall -O3
#FFLAGS      = -O3
FFLAGS += $(OPTIONS)
LDFLAGS     = #-lmpi
LDFLAGS     = #-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/ #-lmpi
endif


ifeq ($(SYSTYPE),"ccage")
F90         = /usr/lib64/openmpi/bin/mpif90#-openmpi-gcc8
#FFLAGS      = -g -fcheck=all -ffree-line-length-none -ffree-form -ftrapv -fbounds-check -fbacktrace
FFLAGS      = -O3 -ffree-line-length-none -ffree-form
FFLAGS += $(OPTIONS)
LDFLAGS     = 
endif

ifeq ($(SYSTYPE),"horizon")
F90         = mpif90
FFLAGS      = -fopenmp -g -traceback -O3 -fpp -ftz  -xavx -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all -heap-arrays10
FFLAGS += $(OPTIONS)
LDFLAGS     =  #-lmpi
endif


#############################################################################
# All objects

# choose the gas mix you want to use
# 1/ HI, D, and dust
#COMPOBJS = module_dust_model.o module_D_1215_model.o module_HI_1216_model.o module_gas_composition_HI_D_dust.o
# 2/ HI and dust
#COMPOBJS = module_dust_model.o module_HI_1216_model.o module_gas_composition_HI_dust.o
# 3/ pure HI
COMPOBJS = module_HI_1216_model.o module_gas_composition_HI.o
# 4/ pure RAMSES hydro variables
#COMPOBJS = module_gas_composition_ramses.o
# 5/ SiII 1190/1193
#COMPOBJS = module_SiII_1190_model.o module_SiII_1193_model.o module_gas_composition_SiII_1190_1193.o
# 6/ SiII 1190/1193 and dust
#COMPOBJS = module_dust_model.o module_SiII_1190_model.o module_SiII_1193_model.o module_gas_composition_SiII_1190_1193_dust.o
# 7/ MgII
#COMPOBJS = module_MgII_2796_model.o module_MgII_2804_model.o module_gas_composition_MgII.o
# 8/ MgII and dust
#COMPOBJS = module_dust_model.o module_MgII_2796_model.o module_MgII_2804_model.o module_gas_composition_MgII_dust.o
# 9/ FeII UV1
#COMPOBJS = module_FeII_2587_model.o module_FeII_2600_model.o module_gas_composition_FeII_UV1.o
# 10/ FeII UV1 and dust
#COMPOBJS = module_dust_model.o module_FeII_2587_model.o module_FeII_2600_model.o module_gas_composition_FeII_UV1_dust.o
# 11/ SiII 1260
#COMPOBJS = module_SiII_1260_model.o module_gas_composition_SiII_1260.o
# 12/ SiII 1260 and dust
#COMPOBJS = module_dust_model.o module_SiII_1260_model.o module_gas_composition_SiII_1260_dust.o
# 13/ pure dust
#COMPOBJS = module_dust_model.o module_gas_composition_dust.o
# 14/ HI, HeI, HeII -> WARNING, this is for fesc measurements only ! 
#COMPOBJS = module_gas_composition_HI_HeI_HeII.o

GASOBJS = module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_ramses.o module_uparallel.o $(COMPOBJS)

OBJSDOM = $(GASOBJS) module_select.o module_mesh.o

OBJSPHO = $(GASOBJS) module_mesh.o module_photon.o 

OBJSMPI = $(OBJSPHO) module_parallel_mpi.o module_worker.o module_master.o

OBJLYAFG = module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_ramses.o module_uparallel.o module_HI_1216_model.o

OBJFESC	= $(GASOBJS) module_mesh.o module_gray_ray.o


#############################################################################
#.SUFFIXES: .o .f90
#.f90.o:
#	$(F90) -c  $(FFLAGS) $<

%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

#############################################################################
all: clean CreateDomDump PhotonsFromSourceModel PhotonsFromStars PhotonsFromHaloes rascas-serial rascas

PhotonsFromSourceModel: $(OBJSPHO) PhotonsFromSourceModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) PhotonsFromSourceModel.o

PhotonsFromStars: $(OBJSPHO) module_ssp_lib.o PhotonsFromStars.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) module_ssp_lib.o PhotonsFromStars.o

PhotonsFromHaloes: $(OBJSPHO) PhotonsFromHaloes.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) PhotonsFromHaloes.o

LyaPhotonsFromGas: $(OBJLYAFG)  LyaPhotonsFromGas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJLYAFG)  LyaPhotonsFromGas.o

CreateDomDump: $(OBJSDOM) CreateDomDump.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) CreateDomDump.o

ExtractSubvol: $(OBJSDOM) ExtractSubvol.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) ExtractSubvol.o

FescFromStars: $(OBJFESC) FescFromStars.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJFESC) FescFromStars.o

rascas-serial: $(OBJSPHO) rascas-serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) rascas-serial.o

rascas: $(OBJSMPI) rascas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI) rascas.o

clean:
	rm -f $(OBJSDOM) $(OBJSPHO) $(OBJSMPI) *.o *.mod
#############################################################################
