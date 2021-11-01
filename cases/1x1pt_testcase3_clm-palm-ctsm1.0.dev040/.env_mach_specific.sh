# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
source /opt/sw/etc/lmod/profile
module -q purge 
module -q load subversion intel-parallel-studio intel netcdf-c netcdf-fortran intel-mpi hdf5 netcdf-c netcdf-fortran
export OMP_STACKSIZE=64M
export MPI_TYPE_DEPTH=16