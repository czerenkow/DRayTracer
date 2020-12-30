cmake_minimum_required(VERSION 3.16)

project(DRayTracer LANGUAGES CXX)

execute_process(COMMAND "date" "--iso-8601=seconds" OUTPUT_VARIABLE BUILD_TIME)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Using versino without ISPC
# TODO: Here there is information how to add ISPC version
# https://github.com/embree/embree/blob/master/tutorials/find_embree/CMakeLists.txt
find_package(embree 3.2 REQUIRED)

# tested with 0.9.9.7
# TODO: Why I can not put explicitly 0.9 or 0.9.9 here? It fails in such case
# TODO: Seems that adding this is not needed. But why? Is it only a header library?
#find_package(glm REQUIRED)

# This is required if exectured with STD C++ paraller algorithm
#find_package(TBB REQUIRED)

# Threads (pthreads?) needed by MPI version   (TODO: refactor it)
# THREADS_PREFER_PTHREAD_FLAG will cause that -pthread rather than -lpthread will be used
# after target_link_libraries Threads::Threads.
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(DRayTracer
    raytracer.cpp
    renderer.cpp
    image.cpp
    random.cpp
    scene.cpp
    main.cpp
)
target_link_libraries(DRayTracer PRIVATE embree)
#target_link_libraries(DRayTracer PRIVATE tbb)
target_link_libraries(DRayTracer PRIVATE OpenMP::OpenMP_CXX)
#target_link_libraries(DRayTracer PRIVATE Threads::Threads)


#################################################################
add_executable(DRayTracerMPI)
target_sources(DRayTracerMPI PRIVATE
    raytracer.cpp
    renderer.cpp
    image.cpp
    random.cpp
    scene.cpp
    main-mpi.cpp
)
target_compile_definitions(DRayTracerMPI PRIVATE BUILD_TIME=${BUILD_TIME})
target_link_libraries(DRayTracerMPI PRIVATE embree)
target_link_libraries(DRayTracerMPI PRIVATE OpenMP::OpenMP_CXX) # In fact we are not using OpenMP in MPI version
target_link_libraries(DRayTracerMPI PRIVATE Threads::Threads)
target_link_libraries(DRayTracerMPI PRIVATE MPI::MPI_CXX)

#################################################################

#add_library(DRayTracerLib STATIC
#  renderer.cpp
#  image.cpp
#  renderer.h
#)
#target_link_libraries(DRayTracerLib PRIVATE embree)
#target_link_libraries(DRayTracerLib PRIVATE OpenMP::OpenMP_CXX)
#target_compile_definitions(DRayTracerLib PRIVATE DRAYTRACERLIB_LIBRARY)

#################################################################

#install(TARGETS DRayTracerMPI
#        CONFIGURATIONS Release
#        RUNTIME DESTINATION /home/universe.dart.spb/pwoloszkiewicz/workspace/docker-mpi/project)
