#----------------------------------------------
# makefile per BGK
# Giorgio Amati, CASPUR 2000-2012
# Giorgio Amati, CINECA 2013-20??
#----------------------------------------------
# notes:
#
# 	1) works correctly using gmake
#	2) check on
#		Intel/AMD
#
#----------------------------------------------
# targets:
#
#		serial ---> serial code
#		openmp ---> serial + openmp
#               mpi    ---> mpi
#		all    ---> mpi + openmp
# 
# flags
#               DBL = -DDOUBLE_P (double precision, default)
#               DBL =            (single precision)
#               * NON BLOCK (or SENDRECV)
#               * DOMAIN DECOMPOSITION (or MULTIBLOCK)
#               * CHANNEL/DUCT/DRIVEN/COUETTE/NOSLIP
#               * FUSED (or ORIGINAL)
#               * PARALLEL (or SERIAL)
#               * FORCING_X (or FORCING_Y/FORCING_Z)
#               * RAW (or HDF5/MPIIO)
#
#
#
#----------------------------------------------
# default  (mpi)
PREPROC = -DPERIODIC -DORIGINAL 
PREPROC = -DCHANNEL -DORIGINAL -DTWO_FIELDS 
PREPROC = -DCHANNEL 
PREPROC = -DORIGINAL -DCOUETTE
PREPROC = -DDUCT 
PREPROC = -DCHANNEL -DORIGINAL 
PREPROC = -DDRIVEN -DORIGINAL
PREPROC = -DDRIVEN 
DBL = -DDOUBLE_P
DBL = 
FIX = 
#TCASE = -DDRIVEN
#OPT = -xhost
#FOPT = -fast -mmic
#COPT = -fast -mmic
#FOPT = -fast 
#COPT = -fast 
#
VER = mpi
LIB =
#
DATE = $(shell date +%d%m%y)
#
# default
#----------------------------------------------
FINE = "********** Fine Compilazione *********"
#----------------------------------------------
#
#
# PGI (default)
SCC = pgcc
SFC = pgf90
CC = mpicc
FC = mpifort
FIX = -DPGI
FOPT = -fast -Mcontiguous -Mnodepchk
COPT = -fast -Mcontiguous -Mnodepchk
OPT  = -fast -Mcontiguous -Mnodepchk
#
ifdef GNU
 FC = gfortran
 CC = gcc
 FIX = -DORIGINAL -DNOMANAGED -DSERIAL 
 FOPT = -O3 -fopenacc
 COPT = -O3
 OPT  = -O3 -fopenacc
 VER  = gnuopenacc
endif
#
ifdef GNU_8
# CC = mpicc
# FC = mpif90
 FC =/m100/home/userinternal/dmolina1/gcc_nvidia/install/bin/gfortran
 CC = gcc
 FIX = -DORIGINAL -DNOMANAGED -DSERIAL -DSERIAL
 FOPT = -O3 -fopenacc
 COPT = -O3
 OPT  = -O3 -fopenacc
 VER  = gnu8
endif
#
ifdef GNU_9
 FC =/m100/home/userinternal/vruggie1/OFFLOAD/install/bin/gfortran
 CC = gcc
 FIX = -DORIGINAL -DNOMANAGED -DSERIAL
 FOPT = -O3 -fopenacc
 COPT = -O3 
 OPT  = -O3 -fopenacc
 VER  = gnu9
endif
#
ifdef GNU_10
 FC =/m100/home/userinternal/vruggie1/OFFLOAD_NEW/install/bin/gfortran
 CC = gcc
 FIX = -DORIGINAL -DNOMANAGED -DSERIAL
 FOPT = -O3 -fopenacc
 COPT = -O3
 OPT  = -O3 -fopenacc
 VER  = gnu10
endif
#
ifdef PWR
 SCC = pgcc 
 SFC = pgf90 -fast -tp=pwr8  -Mcontiguous -Mnodepchk -Minfo=all
 FIX = -DPGI
endif
#
#
ifdef DOUBLE
 DBL = -DDOUBLE_P
endif
#
ifdef COMM
 COMM1 = -DNOBLOCK
else
 COMM1 = -DSENDRECV
endif
#
ifdef IOMPI
 FIX = -DMPIIO
endif
#
#----------------------------------------------
# module
MOD0 = get_mem_mod.o \

MOD1 = precision.o \
      storage.o \
      mod_tools.o \


# 0th level subroutine
OBJ0 = bgk3d.o \

# 1st level subroutine
OBJ1 = setup.o \
       initialize.o \
       diagnostic.o \
       boundaries.o \
       propagation.o \
       save.o\
       finalize.o \

# 2nd level subroutine
OBJ2 = check_case.o \
       input.o \
       alloca.o \
       bcond.obs.o \
       bcond.comm.noblock.packed.o \
       bcond.comm.packed.o \
       bcond.bc.o \
       setup_MPI.o \
       build.bcond.o \
       build.obs.o \
       restore_raw.o \
       init.o\
       hencol.o \
       outdat.o \
       drag.o \
       vtk_3d.binary.o \
       vtk_xy.binary.o \
       vtk_xz.binary.o \
       vtk_yz.binary.o \
       vtk_om.binary.o \
       varm.o \
       prof_k.o \
       prof_j.o \
       prof_i.o \
       flux_X.o \
       flux_Y.o \
       flux_Z.o \
       diagno.o \
       probe.o \
       probe_global.o \
       movef.o\
       save_raw.o\
       profile.o\
       dealloca.o\
       save.mpiio.o \
       restore.mpiio.o  

# 3rd level subroutine
OBJ3 = get_mem.o \
       check_isend.o \
       collision.o \
       col.o \
       col_MC.o \
       read_topology.o


OBJ = $(OBJ1) $(OBJ0) $(OBJ2) $(OBJ3)

#----------------------------------------------
#
all: $(MOD0) $(MOD1) $(OBJ) 
	$(FC) $(OPT) $(MOD1) $(OBJ) $(FIX) $(LIB) -o ../RUN/bgk3d.$(VER).x 
	@echo $(FINE)
#
#----------------------------------------------
#
serial: VER := serial
serial: PREPROC += -DSERIAL 
serial: FC := $(SFC)
serial: CC := $(SCC)
serial: $(MOD0) $(MOD1) $(OBJ) 
	$(FC) $(FOPT) $(DBL) $(MOD1) $(OBJ) $(FIX) -o ../RUN/bgk3d.$(VER).x
	@echo $(FINE)
#
#----------------------------------------------
openacc: FINE = "Under development: works only with pgi compiler!!!"
openacc: FC := pgf90 -acc -Minfo=acc -fast -ta=tesla:cc70,managed -Mcontiguous
openacc: FC := pgf90 -acc -Minfo=acc -fast -ta=tesla:cc70 -Mcontiguous -DNOMANAGED
openacc: CC := pgcc
#openacc: FC := /m100/home/userinternal/dmolina1/gcc_nvidia/install/bin/gfortran -fopenacc
#openacc: CC := /m100/home/userinternal/dmolina1/gcc_nvidia/install/bin/gcc
openacc: PREPROC = -DSERIAL -DPGI -DOPENACC -DDRIVEN 
openacc: VER = openacc
openacc: $(MOD0) $(MOD1) $(OBJ) 
	$(FC) $(OPT) $(DBL) $(MOD1) $(OBJ) -o ../RUN/bgk3d.openacc.x
	@echo $(FINE)
#
#----------------------------------------------
mpi+openacc: FINE = "Under development: works only with pgi compiler!!!"
mpi+openacc: FC := mpif90 -acc -Minfo=acc -fast -ta=tesla:cc70,managed  -Mcontiguous
mpi+openacc: FC := mpif90 -acc -Minfo=acc -fast -ta=tesla:cc70 -Mcontiguous -DNOMANAGE
#mpi+openacc: FC := mpif90 -acc -Minfo=acc -fast -ta=tesla:cc70  -Mcontiguous
mpi+openacc: CC := mpicc
mpi+openacc: PREPROC = -DPGI -DOPENACC -DDRIVEN  -Mnodepchk 
mpi+openacc: VER = mpi+openacc
mpi+openacc: $(MOD0) $(MOD1) $(OBJ) 
	$(FC) $(OPT) $(DBL) $(MOD1) $(OBJ) -o ../RUN/bgk3d.mpi+openacc.x
	@echo $(FINE)

#----------------------------------------------
openmp: OPT += -fopenmp
openmp: FOPT += -fopenmp
openmp: VER := openmp
openmp: serial

#----------------------------------------------
mpi: all
#
#----------------------------------------------
hybrid: FOPT += -mp
hybrid: OPT += -mp
hybrid: VER := hybrid
hybrid: all
#
#----------------------------------------------
%.o %.mod: $(INC) %.f90
	$(FC) $(FOPT) $(DBL) $(PREPROC) $(FIX) $(COMM1) $(TCASE) -c $<

%.o: $(INC) %.F90
	$(FC) $(FOPT) $(DBL) $(PREPROC) $(FIX) $(COMM1) $(TCASE) -c $<

%.o: $(INC) %.f
	$(FC) $(FOPT) $(DBL) $(PREPROC) $(FIX) $(COMM1) $(TCASE) -c $<

%.o: $(INC) %.c
	$(CC) $(COPT)                          -c $<

%.p: %.f90
	$(FC) $(FOPT) $(DBL) $(PREPROC) $(FIX) $(QQQ) $(TCASE) $<

#----------------------------------------------
clean:
	rm -f *.o;
	rm -f *.mod;
	rm -f core;
#----------------------------------------------
info:  $(INC)
	@echo "Modules          =  "$(MOD);
	@echo "Objects          =  "$(OBJ);
	@echo "Compilar         =  "$(FC);
	@echo "Compiler flags   =  "$(OPT);
	@echo "Preproc  flags   =  "$(PREPROC);
	@echo "Fix      flags   =  "$(FIX);
	@echo "Double   flags   =  "$(DBL);
	@echo "Architecture     =  "$(HOSTTYPE);
	@echo "Version          =  "$(VER);
	@echo "Libraries        =  "$(LIB);

#----------------------------------------------
