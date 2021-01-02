# DRayTracer
This is a toy ray tracer that is based on Intel Embree, capable to compute on multiple nodes with MPI.
In the source tree there is also docker-mpi project included (git subtree of https://github.com/czerenkow/docker-mpi) modified in such a way that it installs Embree.

This version is not very flexible as it renders exacly one image :)
![Renderer output](/doc/output.png)


## Build
```
mkdir DRayTracer-build
cd DRayTracer-build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/embree-3.12.1.x86_64.linux /path/to/DRayTracer/
make
# And run on one host:
./DRayTracer
```
Two binaries are built: DRayTracer and DRayTracerMPI. The first is single-host version, and the latter can be run in the MPI environment.

## MPI version
Details how to set up the cluster are here https://github.com/czerenkow/docker-mpi, and this version additionally includes Embree. Only what you need to do is:
```
docker-compose build
docker-compose push
./start-stack
./attach-stack
```
If you need to experiment with Dockerfile, you can download Embree once `cd docker-mpi &&. / Download-embree.sh` and modify Dockerfile (details inside), to use COPY command (to copy downloaded file embree-3.12.1.x86_64 .linux.tar.gz), instead of ADD command (which downloads this file every time).

Now copy the binary `DRayTracerMPI` to `DRayTracer/docker-mpi/project` and then, on docker-mpi claster, you can run the program this way:
```
mpiexec -f machinefile -n 3 project/DRayTracerMPI
```
Output image will be stored in `project` directory.

# Mini benchmark
I have two machines, and the first one is a bit faster than the second one. Both machines are connected by WiFi network (I cannot find my cable ;).

The rendering time of the image image 1944x1920 and 700 samples/pixel:
```
Machine 1: 121 s
Machine 2: 168 s  (slow)
MPI version: 105 s
```

So speedup is respectively 1.15 and 1.6. Far from perfection as perfect speed up should be close to 2. One of the problem in current code is that calculated tiled images are transfered when all calcualtion are done, so network is utilized purely. It will be improved in the next version!



