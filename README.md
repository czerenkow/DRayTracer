# DRayTracer
Toy ray tracer that is besed on Intel Embree and option to be run run on multiple nodes with MPI.
In the source tree there is also docker-mpi project included (git subtree of https://github.com/czerenkow/docker-mpi) modified that installs Embree. This allows to run DRayTracer on multiple nodes in order to speed up calculations.

## Build
```
mkdir DRayTracer-build
cd DRayTracer-build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/embree-3.12.1.x86_64.linux /path/to/DRayTracer/
make
```
Two binaries are built: DRayTracer and DRayTracerMPI. The first is single-host version, and the latter can be run in the MPI environment.

## Setup MPI cluster
Details are described here https://github.com/czerenkow/docker-mpi, however this version includes Embree and this requires manually download this library.
```
cd docker-mpi
./download-embree.sh
```
This downloads embree-3.12.1.x86_64.linux.tar.gz which is required by Dockerfile. All next steps are the same as in not modified docker-mpi. Copy the binary `DRayTracerMPI` to `DRayTracer/docker-mpi/project` and then, on docker-mpi claster, you can run the program this way:
```
mpiexec -f machinefile -n 3 project/DRayTracerMPI
```
