SUPPORTS_CXX := FALSE
ifeq ($(COMPILER),intel)
  FFLAGS_NOOPT :=  -O0 
  SCC :=  icc 
  CXX_LDFLAGS :=  -cxxlib 
  SUPPORTS_CXX := TRUE
  FFLAGS :=  -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source  
  FIXEDFLAGS :=  -fixed  
  CXX_LINKER := FORTRAN
  FC_AUTO_R8 :=  -r8 
  CFLAGS :=   -qno-opt-dynamic-align -fp-model precise -std=gnu99 
  FREEFLAGS :=  -free 
  SFC :=  ifort 
  SCXX :=  icpc 
endif
ifeq ($(COMPILER),gnu)
  FFLAGS_NOOPT :=  -O0 
  SCC :=  gcc 
  HAS_F2008_CONTIGUOUS := FALSE
  SUPPORTS_CXX := TRUE
  FFLAGS :=   -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none 
  FIXEDFLAGS :=   -ffixed-form 
  CXX_LINKER := FORTRAN
  FC_AUTO_R8 :=  -fdefault-real-8 
  CFLAGS :=  -std=gnu99 
  FREEFLAGS :=  -ffree-form 
  SFC :=  gfortran 
  SCXX :=  g++ 
endif
MPIFC :=  mpiifort 
PIO_FILESYSTEM_HINTS := lustre
MPICC :=  mpiicc 
PNETCDF_PATH := 
MPICXX :=  mpiicpc 
NETCDF_PATH := $(NETCDF)
MPI_PATH := $(MPI_ROOT)
MPI_LIB_NAME := mpi
FFLAGS := $(FFLAGS)  -xCORE-AVX2 
CPPDEFS := $(CPPDEFS)  -DCESMCOUPLED 
CPPDEFS := $(CPPDEFS)  -D$(OS) 
SLIBS := $(SLIBS) $(shell nf-config --flibs) -L/usr/users/model4bk/CESM/lapack-3.6.1 -llapack -lrefblas \
     -L/cm/local/apps/gcc/5.2.0/lib64 -lgfortran
ifeq ($(MODEL),micom)
  FFLAGS := $(FFLAGS)  -r8 
endif
ifeq ($(MODEL),fv3gfs)
  FFLAGS := $(FFLAGS)  $(FC_AUTO_R8) 
endif
ifeq ($(MODEL),mom)
  FFLAGS := $(FFLAGS)  $(FC_AUTO_R8) 
endif
ifeq ($(DEBUG),FALSE)
  FFLAGS := $(FFLAGS)  -O2 
endif
ifeq ($(MODEL),fv3gfs)
  CPPDEFS := $(CPPDEFS)  -DSPMD 
endif
ifeq ($(MODEL),pop)
  CPPDEFS := $(CPPDEFS)  -D_USE_FLOW_CONTROL 
endif
ifeq ($(COMPILER),intel)
  CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DCPRINTEL
  ifeq ($(compile_threaded),true)
    FFLAGS := $(FFLAGS)  -qopenmp 
    CFLAGS := $(CFLAGS)  -qopenmp 
  endif
  ifeq ($(DEBUG),TRUE)
    FFLAGS := $(FFLAGS)  -O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created 
    CFLAGS := $(CFLAGS)  -O0 -g 
  endif
  ifeq ($(DEBUG),FALSE)
    FFLAGS := $(FFLAGS)  -O2 -debug minimal 
    CFLAGS := $(CFLAGS)  -O2 -debug minimal 
  endif
  ifeq ($(MPILIB),mvapich2)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpich2)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpt)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),openmpi)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpich)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mvapich)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),impi)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpi-serial)
    SLIBS := $(SLIBS)  -mkl 
  endif
  ifeq ($(compile_threaded),true)
    FFLAGS_NOOPT := $(FFLAGS_NOOPT)  -qopenmp 
    LDFLAGS := $(LDFLAGS)  -qopenmp 
  endif
endif
ifeq ($(COMPILER),gnu)
  CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DNO_R16 -DCPRGNU
  ifeq ($(compile_threaded),true)
    FFLAGS := $(FFLAGS)  -fopenmp 
    CFLAGS := $(CFLAGS)  -fopenmp 
  endif
  ifeq ($(DEBUG),TRUE)
    FFLAGS := $(FFLAGS)  -g -Wall -Og -fbacktrace -ffpe-trap=zero,overflow -fcheck=bounds 
    CFLAGS := $(CFLAGS)  -g -Wall -Og -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds 
  endif
  ifeq ($(DEBUG),FALSE)
    FFLAGS := $(FFLAGS)  -O 
    CFLAGS := $(CFLAGS)  -O 
  endif
  ifeq ($(compile_threaded),true)
    LDFLAGS := $(LDFLAGS)  -fopenmp 
  endif
endif
ifeq ($(MODEL),fv3gfs)
  INCLDIR := $(INCLDIR)  -I$(EXEROOT)/atm/obj/FMS 
endif
