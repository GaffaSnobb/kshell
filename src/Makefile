#
# NOTE:  '-DMPI' option is required for MPI parallel 
# 

# --- Intel Fortran compiler, single node
FC = ifort
FFLAGS = -O2 -qopenmp -no-ipo 
# FFLAGS =  -fast -qopenmp 
# FFLAGS =  -fast -qopenmp -no-ipo -check bounds -check pointers -check uninit -traceback
# FFLAGS =  -g -check all -warn -traceback -debug all -p -no-ipo -qopenmp
# FFLAGS =  -g -check all -warn -traceback -debug all -p -no-ipo
# FFLAGS =  -g -check all -traceback -debug all -p -no-ipo
LIBS = -qmkl # Intel Math Kernel Library for BLAS/LAPACK routines


# --- GNU Fortran compiler, single node 
# FC = gfortran
# FFLAGS = -O3 -fopenmp -fallow-argument-mismatch
# LIBS = -llapack -lblas -lm 



# --- PEGASUS intel mpi ---
# FC = mpiifort -DMPI
# FFLAGS = -O3 -qopenmp
# LIBS = -mkl


# --- PEGASUS nvfortran + openmpi ---
# FC = nvfortran 
# FC = mpif90 -DMPI -DSPARC
# FFLAGS = -mp -O 
# LIBS = -llapack -lblas


# gfortran + openmpi
# FC = mpif90  -DMPI
# FFLAGS = -O3 -fopenmp -fallow-argument-mismatch
# LIBS = -llapack -lblas -lm 
#
# gfortran + mkl
# export MKLROOT=/work/opt/local/cores/intel/compilers_and_libraries_2017.3.191/linux/mkl
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MKLROOT}/lib/intel64
# FC = mpif90 -m64 -I${MKLROOT}/include sample.f90 -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
# 
# MKLROOT=/work/opt/local/cores/intel/compilers_and_libraries_2017.3.191/linux/mkl
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MKLROOT}/lib/intel64
# LIBS = -m64 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#

# Intel MPI
# FC = mpif90 -DMPI -fc=ifort

# --- Fujitsu Fortran compiler, MPI (K computer, FX10, FUGAKU) ---
# FC = mpifrtpx -DMPI
# FFLAGS = -DSPARC -Kopenmp,fast,ocl,ilfunc,preex,noparallel -Nlst=t
# LIBS = -SSL2BLAMP  # K computer (MPI/OpenMP/BLAS,LAPACK(thread version))  recommended

#  FFLAGS =  -DSPARC -Kopenmp,fast,ocl,preex,noparallel -Qt
#  FFLAGS = -DSPARC -Kopenmp,fast,ocl,preex,noparallel,prefetch_indirect,simd=2 \
#     -Qt -Qp -Nrt_tune_func -Nrt_tune_loop
#
#  FFLAGS = -DSPARC -Kopenmp,fast,ocl,visimpact,prefetch_indirect,prefetch_infer,prefetch_cache_level=all,prefetch_strong,prefetch_strong_L2,simd=2,uxsimd -Qt -Qp
#  FFLAGS = -DSPARC -Haefosux -Kopenmp,fast,ocl,ilfunc,preex,noparallel -Qt
#  LIBS = -SSL2 # K computer (single thread)



# --- Fujitsu Fortran compiler single node ---
# FC = frtpx
# FFLAGS = -DSPARC -Kopenmp,fast,ocl,ilfunc,preex,noparallel -Nlst=t
# LIBS = -SSL2BLAMP  


# --- Fujitsu Fortran compiler on Intel (Nagoya CX400)
# FC = mpifrt -DMPI

# --- PGI Fortran compiler ---
#FC = pgf90

# MacOS + ifort 
# FFLAGS =  -fast -openmp -no-ipo -Wl,-stack_size,40000000  # Mac OSX + ifort
# LIBS = -mkl


# FFLAGS := $(FFLAGS) -DNO_POPCNT  # avoid using popcnt, poppar

#USE_ZPARES = 1
USE_ZPARES = 0
ZPARES_DIR = ./zpares


all: kshell.exe transit.exe count_dim.exe
	cp kshell.exe transit.exe count_dim.exe ../bin/

.PHONY: all clean

clean:
	rm -f *.exe *.o *.mod *.lst *__genmod.f90 gmon.out *.a


.SUFFIXES :
.SUFFIXES : .f90 .F90 .f .f95 .F95 .o .mod .c

%.o : %.mod

.f90.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F90.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.f90.mod :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F90.mod :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

class_stopwatch.o: model_space.o 
model_space.o: constant.o
harmonic_oscillator.o: harmonic_oscillator.f90 constant.o rotation_group.o \
  rotation_group.o model_space.o
sp_matrix_element.o: sp_matrix_element.f90 harmonic_oscillator.o \
  rotation_group.o constant.o
interaction.o: rotation_group.o model_space.o harmonic_oscillator.o sp_matrix_element.o \
  operator_jscheme.o operator_mscheme.o
operator_jscheme.o: model_space.o harmonic_oscillator.o
operator_mscheme.o: model_space.o operator_jscheme.o rotation_group.o
lanczos.o: constant.o lib_matrix.o model_space.o  operator_mscheme.o \
  bridge_partitions.o bp_io.o
partition.o: model_space.o lib_matrix.o class_stopwatch.o constant.o
wavefunction.o: constant.o model_space.o partition.o class_stopwatch.o
bridge_partitions.o: partition.o wavefunction.o \
   operator_mscheme.o model_space.o constant.o class_stopwatch.o
bp_io.o: bridge_partitions.o interaction.o
bp_expc_val.o: bridge_partitions.o
kshell.o: bridge_partitions.o lanczos.o class_stopwatch.o rotation_group.o \
 interaction.o bp_io.o bp_expc_val.o bp_block.o  block_lanczos.o
transit.o: model_space.o interaction.o operator_mscheme.o partition.o \
 bridge_partitions.o class_stopwatch.o rotation_group.o bp_io.o bp_expc_val.o
transit_em.o: model_space.o interaction.o operator_mscheme.o partition.o \
 bridge_partitions.o class_stopwatch.o rotation_group.o bp_io.o  bp_expc_val.o
mup_operate.o: model_space.o interaction.o operator_mscheme.o partition.o \
 bridge_partitions.o class_stopwatch.o rotation_group.o bp_io.o


bp_block.o:  bp_block.F90 bp_block_inc.F90 bridge_partitions.o
block_lanczos.o: wavefunction.o bp_block.o bp_io.o lanczos.o

convert_m2j.exe: model_space.o interaction.o operator_mscheme.o partition.o \
	class_stopwatch.o rotation_group.o bp_io.o

print_jwav.o: model_space.o partition.o rotation_group.o


OBJS = model_space.o interaction.o harmonic_oscillator.o constant.o \
	rotation_group.o sp_matrix_element.o \
	operator_jscheme.o operator_mscheme.o lib_matrix.o lanczos.o \
	partition.o  wavefunction.o  bridge_partitions.o \
	bp_io.o bp_expc_val.o class_stopwatch.o bp_block.o block_lanczos.o


ifdef ZPARES
 FFLAGS := $(FFLAGS) -DZPARES -I$(ZPARES_DIR)/include 
 LIBS := -L$(ZPARES_DIR)/lib  -lzpares  $(LIBS)
 OBJS := $(OBJS)  ss_method.o  eigdensity_estimator.o
 kshell.o: ss_method.o
 ss_method.o: wavefunction.o bp_block.o eigdensity_estimator.o block_lanczos.o
 eigdensity_estimator.o:  class_stopwatch.o
endif



kshell.exe: kshell.o $(OBJS)
	$(FC) $(FFLAGS) -o $@ kshell.o $(OBJS) $(LIBS) 

transit.exe: transit.o $(OBJS)
	$(FC) $(FFLAGS) -o $@ transit.o $(OBJS) $(LIBS) 

transit_em.exe: transit_em.o $(OBJS)
	$(FC) $(FFLAGS) -o $@ transit_em.o $(OBJS) $(LIBS) 

mup_operate.exe: mup_operate.o $(OBJS)
	$(FC) $(FFLAGS) -o $@ mup_operate.o $(OBJS) $(LIBS) 

count_dim.exe: count_dim.f90 $(OBJS)
	$(FC) $(FFLAGS)  -static -o $@ count_dim.f90 $(OBJS) $(LIBS) 


overlap_rwf.exe: overlap_rwf.F90 input_nml.o $(OBJS)
	$(FC) $(FFLAGS) -o $@ overlap_rwf.F90 input_nml.o $(OBJS) $(LIBS)

print_mbit_wf.exe: print_mbit_wf.f90 $(OBJS)
	$(FC) $(FFLAGS) -o $@ print_mbit_wf.f90 $(OBJS) $(LIBS)



convert_m2j.exe: convert_m2j.F90 $(OBJS)
	$(FC) $(FFLAGS) -o $@ convert_m2j.F90 $(OBJS) $(LIBS) 


print_jwav.exe: print_jwav.o $(OBJS)
	$(FC) $(FFLAGS) -o $@ print_jwav.F90 $(OBJS) $(LIBS) \
	&& cp print_jwav.exe ../bin/

