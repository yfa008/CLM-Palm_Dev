SUPPORTS_CXX := FALSE
ifeq ($(COMPILER),intel)
  FFLAGS :=  -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source  
  SUPPORTS_CXX := TRUE
  CFLAGS :=   -qno-opt-dynamic-align -fp-model precise -std=gnu99 
  CXX_LDFLAGS :=  -cxxlib 
  CXX_LINKER := FORTRAN
  FC_AUTO_R8 :=  -r8 
  FFLAGS_NOOPT :=  -O0 
  FIXEDFLAGS :=  -fixed  
  FREEFLAGS :=  -free 
endif
MPICC :=  mpiicc  
MPICXX :=  mpiicpc 
MPIFC :=  mpiifort 
SCC :=  icc 
SCXX :=  icpc 
SFC :=  ifort 
SLIBS :=  -L$(NETCDF_PATH)/lib -lnetcdf -lnetcdff  -lhdf5_hl -lhdf5 -ldl -lm -lz -lcurl -mkl=cluster  
NETCDF_PATH :=  $(NETCDF_PATH)
CPPDEFS := $(CPPDEFS)  -DCESMCOUPLED 
ifeq ($(MODEL),pop)
  CPPDEFS := $(CPPDEFS)  -D_USE_FLOW_CONTROL 
endif
ifeq ($(MODEL),fv3gfs)
  CPPDEFS := $(CPPDEFS)  -DSPMD 
  FFLAGS := $(FFLAGS)  $(FC_AUTO_R8) 
endif
ifeq ($(MODEL),mom)
  FFLAGS := $(FFLAGS)  $(FC_AUTO_R8) 
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
  ifeq ($(MPILIB),mpich)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpich2)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mvapich)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mvapich2)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpt)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),openmpi)
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
ifeq ($(MODEL),fv3gfs)
  INCLDIR := $(INCLDIR)  -I$(EXEROOT)/atm/obj/FMS 
endif
