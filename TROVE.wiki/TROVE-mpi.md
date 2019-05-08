# MPI functionality

Distributed memory parallelism using MPI has been implemented in the
PTcontracted_matalem_class subroutine, plus subroutines called by or related to
it. Most MPI functionality has been consolidated in the mpi_aux.f90 file.

File output has been implemented using MPI-IO. Because MPI is not compatible
with the Fortran's record-based unformatted file format, a new file format has
been implemented. Preserving the old file format by using Master I/O has been
considered, but cannot be supported for large simulations (where distributed
memory parallelism is required) as the size of the data needing to be written
will exceed the available memory on a compute node.

The new file format is similar to the non-mpi format, except for using 'raw'
binary data versus Fortran's records. Conversion utilities between the two
formats will be provided. 

The following has been implemented:

 - Writing of single-file data
 - Reading of single-file data
 - Writing of split file data

The following needs testing:

 - Reading of split file data

The following is in progress:

 - Reading the new format on MPI-free systems
 - Conversion functionality


To enable writing in the new file format, add the parameter `format		mpiio` to
the input file under `CHECK_POINT`, e.g.:

```
CHECK_POINT
HAMILTONIAN read
potential   read
kinetic     read
external    read
basis_set   read
CONTRACT    read
matelem     save
extmatelem  save
eigenfunc   save
format			mpiio
END
```

NOTE: As other routines than `matelem` have not been touched yet, expect a large
amount of duplicated output from stdout. This does not affect the integrity of
file output.

# How to run

## Requirements

The only additional requirement to build and run the MPI version of TROVE is a
compatible MPI library. The code has been tested with Intel MPI, versions
2017.4, 2018.4 and 2019.3. Intel 2018.1 is confirmed *not* working.

This version of TROVE should also build and run correctly with OpenMPI, and
should be compatible with GNU GFortran but this is untested.

Although not required, it is highly recommended to run on a distributed
filesystem (e.g. Lustre + striping) for significantly increased I/O performance.

## Building TROVE-MPI

Make sure the compiler is set to 'mpif90' or 'mpifort' in the Makegfile, then
run `make` in the source directory.

## Running TROVE-MPI

Once built, run the code like the sequential version but prepended with `mpirun`
or your scheduler's mpi wrapper (e.g. `srun`). E.g.:

`mpirun trove.x infile.inp`

A sample batch script for SLURM is provided.
