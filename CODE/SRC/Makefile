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
#               FERMI: FIX = -WF,-DDEBUG_1,-DDEBUG_2,-DFERMI_FIX
#               KNL: make hybrid "FIX = -qopt-report-phase=vec -xmic-avx512 -qopt-prefetch=5 -qopt-report=4 -DKNOWN_TRIP_COUNT"
#               KNL: make openmp KNL=1  " FIX=-assume contiguous_pointer -qopt-zmm-usage=high " 
#                Pointers --> FIX=-assume contiguous_pointer 
#
#
#
#----------------------------------------------
# default  (mpi)
PREPROC = -DDRIVEN 
DBL = 
FIX = 
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
# generic option
FOPT = -O2
COPT = -O2
OPT = -O2
#
# generic compiler
CC = mpiicc
FC = mpiifort
SCC = icc
SFC = ifort
#
ifdef PGI
 SCC = pgcc
 SFC = pgf90
 CC = mpicc
 FC = mpifort
 FC = mpipgifort
 FIX = -DPGI
 FOPT = -fast -Mcontiguous -Mnodepchk
 COPT = -fast -Mcontiguous -Mnodepchk
 OPT = -fast -Mcontiguous -Mnodepchk
endif
#
ifdef FUGAKU
 SFC = frtpx
 SCC = FCCpx
 FC = mpifrtpx
 CC = mpiFCCpx
endif
#
ifdef GNU
 CC = mpicc
 FC = mpif90
 SFC = gfortran
 SCC = gcc
endif
#
ifdef ARM
 CC = mpicc
 FC = mpif90
 SFC = armflang
 SCC = armclang
# CC = gcc
# FC = gfortran
# SFC = gfortran
# SCC = gcc
# FOPT = -O3 -mcpu=native  -march=armv8-a+sve 
 FOPT = -Ofast -fno-stack-arrays -mcpu=native
 LIB = -L/opt/arm/gcc-8.2.0_Generic-AArch64_RHEL-7_aarch64-linux/lib64/
# FOPT = -Ofast -march=armv8-a+sve
endif
#
ifdef BDW
 FOPT = -O3 -xCORE-AVX2 -assume contiguous_pointer
 COPT = -O3 -xCORE-AVX2
endif
#
ifdef KNL
 FOPT = -O3 -xMIC-AVX512 -qopt-prefetch=5 -assume contiguous_pointer
 COPT = -O3 -xMIC-AVX512 -qopt-prefetch=5
endif
#
ifdef SKL
 FOPT = -O3 -xCORE-AVX512 -mtune=skylake-avx512 -assume contiguous_pointer
 COPT = -O3 -xCORE-AVX512 -mtune=skylake-avx512
# FOPT = -O3 -mtune=skylake-avx512
# COPT = -O3 -mtune=skylake-avx512
endif
#
ifdef PWR
 SCC = pgcc 
 SFC = pgf90 -fast -tp=pwr8  -Mcontiguous -Mnodepchk -Minfo=all
 FIX = -DPGI
endif
#

ifdef LIKWID
 FIX = -DLIKWID -I/cineca/tmp/likwid_demON/include/ -L/cineca/tmp/likwid_demON/lib/ -llikwid
 FINE = " using LIKWID API  only marconi"
endif

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
ifdef IOHDF5
 FIX =
endif
# include file
INC = 
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
       bcond.comm.o \
       bcond.comm.openmp.o \
       bcond.comm.noblock.o \
       bcond.comm.noblock.packed.o \
       bcond.comm.noblock.packed.try.o \
       bcond.comm.noblock.packed.overlap.o \
       bcond.comm.packed.o \
       bcond.bc.o \
       bcond.driven.o \
       setup_MPI.o \
       build.bcond.o \
       build.obs.o \
       restore_raw.o \
       init.o\
       hencol.o \
       outdat.o \
       drag.o \
       vtk_3d.o \
       vtk_3d.binary.o \
       vtk_xy.o \
       vtk_xy.binary.o \
       vtk_xz.o \
       vtk_xz.binary.o \
       vtk_yz.o \
       vtk_yz.binary.o \
       vtk_om.o \
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
       wall.o \
       slab.o \
       sphere.o \
       movef.o\
       movef_TF.o\
       save_raw.o\
       profile.o\
       dealloca.o\
       save.mpiio.o \
       restore.mpiio.o  
#       bcond.comm.fullduplex.o  \

# 3rd level subroutine
OBJ3 = get_mem.o \
       check_isend.o \
       collision.o \
       col.o \
       col_MC.o \
       col_MC.box.o \
       read_topology.o

# openacc version
#OBJ5 = all_in_one.o\

# hdf5 version
OBJ6 = save.hdf5.parallel.o \
       save.hdf5.serial.o  \
       save.obs.hdf5.serial.o  \
       hdf5_2d_xz.serial.o \
       hdf5_2d_xz.parallel.o \
       hdf5_2d_xy.parallel.o \
       restore.hdf5.serial.o \
       restore.obs.hdf5.parallel.o \
       restore.obs.hdf5.serial.o \
       restore.hdf5.parallel.o 
   
# cudafortran version
OBJ7 = module_CUDA.o \
       collision.o \
       col.o 

#OBJ7 = all_in_one.cudaf.o \

OBJ = $(OBJ1) $(OBJ0) $(OBJ2) $(OBJ3)

#----------------------------------------------
#
all: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(FIX) $(LIB) -o ../RUN/bgk3d.$(VER).x 
	@echo $(FINE)

#----------------------------------------------
# my laptop
hdf5: FC := /home/giorgio/LOCAL_SW/hdf5_serial/bin/h5fc -DSERIAL
hdf5: FC := h5pfc -qopenmp
hdf5: CC := mpiicc
hdf5: PREPROC += -DHDF5 
#hdf5: FIX := -DDEBUG_1
#hdf5: DBL := 
#hdf5: FC := h5fc
#hdf5: CC := icc
#hdf5: PREPROC += -DSERIAL 
hdf5: VER := hdf5
hdf5: $(MOD0) $(MOD1) $(OBJ) $(OBJ4) $(OBJ6)
	$(FC) $(OPT) $(DBL) $(MOD1) $(OBJ) $(OBJ4) $(OBJ6) -o ../RUN/bgk3d.$(VER).x
	@echo $(FINE)

#----------------------------------------------
try: FINE = "Under development: check preprocessing flags!!!"
try: FC := mpif90
try: CC := mpiicc
#try: FC := gfortran
#try: CC := gcc
try: PREPROC = -DORIGINAL -DREGOLARIZED -DCHANNEL
try: PREPROC = -DORIGINAL -DREGOLARIZED -DDRIVEN -DSERIAL
#try: PREPROC = -DORIGINAL -DPERIODIC
#try: DBL = 
try: VER := try
try: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(DBL) $(MOD1) $(OBJ) $(OBJ4) -o ../RUN/bgk3d.$(VER).x
	@echo $(FINE)

#----------------------------------------------
mb: FINE = "Multi block development: check preprocessing flags!!!"
mb: FC := mpif90
mb: CC := mpicc
mb: PREPROC += -DMB
mb: VER := mb
mb: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(DBL) $(MOD1) $(OBJ) $(OBJ4) -o ../RUN/bgk3d.$(VER).x
	@echo $(FINE)


#----------------------------------------------
#serial: FC := pgf90 -Minfo=all 
#serial: CC := pgcc
#serial: PREPROC += -DSERIAL -DPGI
serial: VER := serial
serial: PREPROC += -DSERIAL 
serial: FC := $(SFC)
serial: CC := $(SCC)
serial: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(FOPT) $(DBL) $(MOD1) $(OBJ) $(OBJ4) $(FIX) $(TCASE)-L/opt/arm/gcc-8.2.0_Generic-AArch64_RHEL-7_aarch64-linux/lib64/ -o ../RUN/bgk3d.$(VER).x
	@echo $(FINE)

#----------------------------------------------
couette: FINE = "Couette Flow: under development!!!"
couette: FC := mpif90 -O3
couette: CC := mpiicc
couette: PREPROC = -DCOUETTE -DORIGINAL 
couette: DBL :=
couette: TCASE = -DCOUETTE
couette: VER = couette
couette: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(DBL) $(MOD1) $(OBJ) $(OBJ4) -o ../RUN/bgk3d.$(VER).x
	@echo $(FINE)


#----------------------------------------------
openacc: FINE = "Under development: works only with pgi compiler!!!"
openacc: FC := pgf90 -fast -ta=multicore -Mcontiguous
openacc: FC := pgf90 -acc -Minfo=acc -fast -ta=tesla:cc70 -Mcontiguous
openacc: FC := pgf90 -acc -Minfo=acc -fast -ta=tesla:cc35,managed -Mcontiguous
openacc: FC := pgf90 -acc -Minfo=acc -fast -ta=tesla:cc70,managed -Mcontiguous
openacc: CC := pgcc
#openacc: PREPROC = -DSERIAL -DORIGINAL -DPGI -DOPENACC -DDRIVEN 
openacc: PREPROC = -DSERIAL -DPGI -DOPENACC -DDRIVEN 
#openacc: PREPROC = -DSERIAL -DPGI -DOPENACC -DDUCT  -Mnodepchk
openacc: VER = openacc
openacc: $(MOD0) $(MOD1) $(OBJ) $(OBJ5)
	$(FC) $(OPT) $(DBL) $(MOD1) $(OBJ) $(OBJ5) -o ../RUN/bgk3d.openacc.x
	@echo $(FINE)

#----------------------------------------------
cudaf: FINE = "Under development: check preprocessing flags!!!"
cudaf: FC := pgf90 -Mcuda=5.0,cc3.5 -Minfo=all
#cudaf: FC := pgf90 -Mcuda=cc2.0 -Minfo=all
cudaf: CC := pgcc
cudaf: PREPROC = -DSERIAL -DORIGINAL -DPGI -DCUDAFOR -DDRIVEN
cudaf: PREPROC = -DSERIAL -DPGI -DCUDAFOR -DPERIODIC
#cudaf: DBL := 
cudaf: VER = cudaf
cudaf: precision.o module_CUDA.o $(MOD0) $(MOD1) $(OBJ) $(OBJ7)
	$(FC) $(OPT) $(DBL) $(MOD1) $(OBJ) $(OBJ7) -o ../RUN/bgk3d.cudaf.x
	@echo $(FINE)



# "FC=pgf90 -Mcuda -Minfo" "CC=pgcc" "FIX=-DSERIAL -DPGI -DCUDAFOR -DORIGINAL" "DBL = " 

#----------------------------------------------
mpi+openacc: FINE = "Under development: works only with pgi compiler!!!"
#mpi+openacc: FC := mpif90 -acc -Minfo=acc -fast -ta=tesla:cc60,managed  -Mcontiguous
mpi+openacc: FC := mpif90 -acc -Minfo=acc -fast -ta=tesla:cc70,managed  -Mcontiguous
#mpi+openacc: FC := mpipgifort -acc -Minfo=acc -fast -ta=tesla:cc70,managed  -Mcontiguous
mpi+openacc: CC := mpicc
#mpi+openacc: CC := mpipgicc
mpi+openacc: PREPROC = -DPGI -DOPENACC -DDRIVEN  -Mnodepchk -DCOMM=1
mpi+openacc: VER = mpi+openacc
mpi+openacc: $(MOD0) $(MOD1) $(OBJ) $(OBJ5)
	$(FC) $(OPT) $(DBL) $(MOD1) $(OBJ) $(OBJ5) -o ../RUN/bgk3d.mpi+openacc.x
	@echo $(FINE)

#----------------------------------------------
#openmp: OPT +=  -qopenmp
#openmp: FOPT += -qopenmp
#openmp: FOPT += -mp
#openmp: OPT += -mp
openmp: OPT += -fopenmp
openmp: FOPT += -fopenmp
openmp: VER := openmp
openmp: serial

#----------------------------------------------
mpi: all

#----------------------------------------------
#hybrid: FOPT += -fopenmp
#hybrid: OPT += -fopenmp
#hybrid: FOPT += -mp
#hybrid: OPT += -mp
hybrid: FOPT += -mp
hybrid: OPT += -mp
hybrid: VER := hybrid
hybrid: all

#----------------------------------------------
mic: FC := mpifc -mmic 
mic: CC := mpiicc -mmic 
mic: FC := mpifc -mmic -openmp
mic: CC := mpiicc -mmic -openmp
mic: FOPT = -fast
mic: OPT = -fast
mic: VER := mic
mic: all

#----------------------------------------------
#mpip: LIB += -L../UTIL/MPIP/mpiP-3.4.1-gnu/lib/  -lmpiP -lm -lbfd -liberty
mpip: LIB += -L../UTIL/MPIP/mpiP-3.4.1-intel/lib/  -lmpiP -lm -lbfd -liberty
mpip: VER := mpiP
mpip: all

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
prep: FIX := -E
prep: col.p
	@echo "only preprocessing "
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
test:
	make clean;
	make mpi;
	make clean;
	make serial
	make clean;
	make hdf5;
	make clean;
	make openmp;
	make clean;
	make hybrid;
	make clean;
	make mic;
	make clean;
	make openacc;
#----------------------------------------------
tar:
	make clean
	rm -rf ../RUN/*.vtk
	rm -rf ../RUN/*.bin
	rm -rf ../RUN/*.h5
#	cd ../TEST/
#	./clean.x
#	cd -
	tar -cvf ../../archive.bgk3d.$(DATE).tar ../../BGK_3D/*
#----------------------------------------------
help:
	@echo "different target available"
	@echo "mpi (default)"
	@echo "serial (raw I/O)"
	@echo "hdf5  (paralle+hdf5)"
	@echo "openmp (serial+openmp)"
	@echo "hybrid (mpi+openmp)"
	@echo "openacc (serial+openac)"
	@echo "mic (mpi for intel MIC)"
	@echo "couette"
	@echo "try (regularized, not working"
	@echo "cudaf"
	@echo "mb (multiblock, not working)"
#----------------------------------------------

