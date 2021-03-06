#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <sys/sysinfo.h> // number of cpu

//#include <execution>
//#include <tbb/tbb.h>

#include "renderer.h" // inteface to this module

#include "image.h"
#include "utils.h"
#include "raytracer.h"
#include "random.h"
#include "constants.h"
#include "camera.h"
#include "scene.h"
#include "thread-pool.h"

using std::cout;
using std::endl;




struct Renderer {
    Image image{frame_columns, frame_rows};
    CameraRegular camera{{-3,7,2}, {0,0,0}, {0,0,1}, 40, float{frame_columns} / frame_rows /* aspect */};

    RTCDevice device;
    std::unique_ptr<RayTracer> rt;


    glm::vec3 getSample(float xf_rand, float yf_rand)
    {
        RTCRayHit rayh = initRayHit();
        camera.getRay(xf_rand/frame_columns, yf_rand/frame_rows, rayh.ray);
        return rt->rayColor(rayh);
    }

    void renderRegion(std::size_t x_start, std::size_t x_stop,
                      std::size_t y_start, std::size_t y_stop)
    {
        for (std::size_t yn = y_start; yn < y_stop; yn++) {
            for (std::size_t xn = x_start; xn < x_stop; xn++) {
                const float xf = static_cast<const float>(xn);
                const float yf = static_cast<const float>(yn);

                glm::vec3 color {0.0};
                for (int s = 0; s < samples; s++) {
                    const float yf_rand = yf + randomns::get01();
                    const float xf_rand = xf + randomns::get01();
                    color += getSample(xf_rand, yf_rand);
                }
                color /= samples;

                color = glm::vec3{glm::clamp(0.0f, 1.0f, color.x),
                        glm::clamp(0.0f, 1.0f, color.y),
                        glm::clamp(0.0f, 1.0f, color.z)};
                color = glm::vec3{std::sqrt(color.x), std::sqrt(color.y), std::sqrt(color.z)};

                // Sample position (0,0) is top left corner of mage, but image.setColor(0,0) is
                // position of top left corner. We need to invert Y axis here.
                image.setColor(xn, frame_rows - 1 - yn, color);
            }
        }
    }

    void renderTile(std::size_t tile) {
        TileInfo t{tile};
        renderRegion(t.c_start, t.c_stop, t.r_start, t.r_stop);
    }

    void renderAllWithOMP() {
        #pragma omp parallel for schedule(dynamic, 1)
        for (std::size_t tile = 0; tile < tiles_rows*tiles_columns; tile++) {
            renderTile(tile);
        }
        image.writeToFileBMP("output.bmp");
    }

    void renderAllWithThreadPool() {
        ThreadPool pool{get_nprocs()};
        for (std::size_t tile = 0; tile < tiles_rows*tiles_columns; tile++) {
            pool.pushTask([this, tile](int){
                renderTile(tile);
            });
        }
        pool.shutdown();

        // Different experiments with parallelisation
        //    std::vector<std::size_t> tiles_list(tiles_rows*tiles_columns);
        //    for (std::size_t i = 0; i < tiles_list.size(); i++) {
        //        tiles_list[i] = i;
        //    }
        //tbb::parallel_for(std::size_t(0), std::size_t(tiles_rows*tiles_columns), )
        //std::for_each(std::execution::par, std::begin(tiles_list),  std::end(tiles_list), render_tile);
        image.writeToFileBMP("output.bmp");
    }

    void writeImage(const char* path) {
        image.writeToFileBMP(path);
    }

    Renderer() {
        device = rtcNewDevice("threads=1,isa=avx");
        if (!device) {
            throw "rtcNewDevice";
        }
        //RTCScene scene = cubeX2InstanceScene(device);
        //RTCScene scene = cubeScene(device);
        //RTCScene scene = spheresScene(device);
        auto [scene, geometry_id, data] = spheresSmallScene(device);
        rt = std::make_unique<RayTracer>(scene, std::move(data));
    }

    ~Renderer() {
        // TODO: free buffers?
        // Release device and all referenced things.
        rtcReleaseDevice(device);
    }
};

//========================================
// class DRenderer
// Simple interface to Renderer class.
//========================================

DRenderer::DRenderer(): r{new Renderer{}}
{
}

DRenderer::~DRenderer() {
    delete r;
}

void DRenderer::renderAllTiles() {
    r->renderAllWithOMP();
    //r->renderAllWithThreadPool();
}

void DRenderer::renderTile(int tile_id) {
    r->renderTile(tile_id);
}

void DRenderer::writeImage() {
    r->writeImage("output.bmp");
}

Image& DRenderer::getImageRef() {
    return r->image;
}

