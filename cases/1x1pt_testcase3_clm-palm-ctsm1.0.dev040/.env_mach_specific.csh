# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
None purge 
None load subversion intel-parallel-studio intel netcdf-c netcdf-fortran intel-mpi hdf5 netcdf-c netcdf-fortran
setenv OMP_STACKSIZE 64M
setenv MPI_TYPE_DEPTH 16