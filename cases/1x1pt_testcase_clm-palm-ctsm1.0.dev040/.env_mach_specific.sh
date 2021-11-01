source /cm/local/apps/environment-modules/3.2.10/init/sh
module purge 
module load subversion intel/compiler intel/mkl netcdf/intel/64 intel/mpi hdf5/parallel/intelmpi/1.8.11 netcdf/intel/parallel
export OMP_STACKSIZE=64M
export MPI_TYPE_DEPTH=16