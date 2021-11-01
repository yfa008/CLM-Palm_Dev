# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
source /cm/local/apps/environment-modules/3.2.10/init/csh
module purge 
module load subversion intel/compiler intel/mkl netcdf/intel/64 intel/mpi hdf5/parallel/intelmpi/1.8.11 netcdf/intel/parallel
setenv OMP_STACKSIZE 64M
setenv MPI_TYPE_DEPTH 16