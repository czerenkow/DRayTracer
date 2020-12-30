#include <iostream>
#include <cstdlib>

#include <limits>
#include <tuple>
#include <chrono>
#include <random>
#include <algorithm>
#include <memory>

//#include <execution>
//#include <tbb/tbb.h>

#include "image.h"
#include "utils.h"
#include "raytracer.h"
#include "random.h"
#include "constants.h"

using std::cout;
using std::endl;


std::ostream& operator<<(std::ostream& os, const glm::vec3& v) {
    os << '[' << v.x << ',' << v.y << ',' << v.z << ']';
    return os;
}




/// \brief The Camera class
/// Position of the screen is on plane z=-1 in camera space (basis). Rays goes from origin through the screen.
class CameraRegular {
public:
    /// \arg vfov is top to bottom in degrees
    /// \arg vup vector up
    /// \arg aspect width/height
    CameraRegular(const glm::vec3& look_from,
                  const glm::vec3& look_at,
                  const glm::vec3& vup,
                  const float vfov,
                  const float aspect):
        origin{look_from}
    {
        // Ortonormal basis of camera
        const auto w = glm::normalize(look_from - look_at); // --> z
        const auto u = glm::normalize(glm::cross(vup, w));  // --> x
        const auto v = glm::cross(w,u);                     // --> y

        // Screen size
        const float theta = vfov * glm::pi<float>() / 180;
        const float half_height = glm::tan(theta / 2);
        const float half_width = aspect * half_height;

        // Screen in 3D in terms of vectors
        lower_left_corner_point = origin - half_width*u - half_height*v - w;
        horizontal = 2 * half_width * u;
        vertical = 2 * half_height * v;
    }

    /// 0 <= x,y <= 1
    /// (0,0) --> lower left corner
    /// (1,1) --> top right corner
    void getRay(float x, float y, RTCRay& ray) const {
        ray.org_x = origin.x;
        ray.org_y = origin.y;
        ray.org_z = origin.z;
        const auto direction = lower_left_corner_point + x*horizontal + y*vertical - origin;
        ray.dir_x = direction.x;
        ray.dir_y = direction.y;
        ray.dir_z = direction.z;
    }

    const glm::vec3 origin;
    glm::vec3 lower_left_corner_point;
    glm::vec3 horizontal; // lenght == horizontal screen size
    glm::vec3 vertical;   // lenght == vertical screen size
};


struct Data {
  glm::vec4* face_colors;
  glm::vec4* vertex_colors;

  Data(std::size_t vertex_n, std::size_t face_n) {
      vertex_colors = static_cast<glm::vec4*>(aligned_alloc(16, vertex_n*sizeof(glm::vec4)));
      face_colors = static_cast<glm::vec4*>(aligned_alloc(16, face_n*sizeof(glm::vec4)));
  }

  ~Data() {
      if (vertex_colors) free(vertex_colors);
      if (face_colors) free(face_colors);
  }

  Data(Data&& data) {
      face_colors = data.face_colors;
      vertex_colors = data.vertex_colors;
      data.face_colors = nullptr;
      data.vertex_colors = nullptr;
  }

  Data& operator=(const Data&) = delete;
  Data(Data& data) = delete;
};



struct DataVertex {
    glm::vec4* vertex;

    DataVertex(std::size_t vertex_n) {
        vertex = static_cast<glm::vec4*>(aligned_alloc(16, vertex_n*sizeof(glm::vec4)));
    }

    ~DataVertex() {
        if (vertex) free(vertex);
    }

    DataVertex(DataVertex&& data) {
        vertex = data.vertex;
        data.vertex = nullptr;
    }

    DataVertex& operator=(DataVertex&& data) {
        vertex = data.vertex;
        data.vertex = nullptr;
        return *this;
    }


    DataVertex& operator=(const DataVertex&) = delete;
    DataVertex(DataVertex& data) = delete;
};





RTCGeometry cubeGeometry(RTCDevice device)
{
  /* create a triangulated cube with 12 triangles and 8 vertices */
  RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

  /* set vertices and vertex colors */
  glm::vec4* vertices = (glm::vec4*) rtcSetNewGeometryBuffer(mesh,
                                                             RTC_BUFFER_TYPE_VERTEX,
                                                             0 /* slot */,
                                                             RTC_FORMAT_FLOAT3,
                                                             sizeof(glm::vec4) /* byteStride */,
                                                             8 /* itemCount */);
  Data data{8 /*vertices*/, 12 /*faces*/};
  data.vertex_colors[0] = glm::vec4(0,0,0, 0); vertices[0].x = -1; vertices[0].y = -1; vertices[0].z = -1;
  data.vertex_colors[1] = glm::vec4(0,0,1, 0); vertices[1].x = -1; vertices[1].y = -1; vertices[1].z = +1;
  data.vertex_colors[2] = glm::vec4(0,1,0, 0); vertices[2].x = -1; vertices[2].y = +1; vertices[2].z = -1;
  data.vertex_colors[3] = glm::vec4(0,1,1, 0); vertices[3].x = -1; vertices[3].y = +1; vertices[3].z = +1;
  data.vertex_colors[4] = glm::vec4(1,0,0, 0); vertices[4].x = +1; vertices[4].y = -1; vertices[4].z = -1;
  data.vertex_colors[5] = glm::vec4(1,0,1, 0); vertices[5].x = +1; vertices[5].y = -1; vertices[5].z = +1;
  data.vertex_colors[6] = glm::vec4(1,1,0, 0); vertices[6].x = +1; vertices[6].y = +1; vertices[6].z = -1;
  data.vertex_colors[7] = glm::vec4(1,1,1, 0); vertices[7].x = +1; vertices[7].y = +1; vertices[7].z = +1;

  /* set triangles and face colors */
  int tri = 0;
  glm::uvec4* triangles = (glm::uvec4*) rtcSetNewGeometryBuffer(mesh,
                                                                RTC_BUFFER_TYPE_INDEX,
                                                                0,
                                                                RTC_FORMAT_UINT3,
                                                                sizeof(glm::uvec4),
                                                                12);
  // left side
  data.face_colors[tri] = glm::vec4(1,0,0, 0); triangles[tri].x = 0; triangles[tri].y = 1; triangles[tri].z = 2; tri++;
  data.face_colors[tri] = glm::vec4(1,0,0, 0); triangles[tri].x = 1; triangles[tri].y = 3; triangles[tri].z = 2; tri++;

  // right side
  data.face_colors[tri] = glm::vec4(0,1,0, 0); triangles[tri].x = 4; triangles[tri].y = 6; triangles[tri].z = 5; tri++;
  data.face_colors[tri] = glm::vec4(0,1,0, 0); triangles[tri].x = 5; triangles[tri].y = 6; triangles[tri].z = 7; tri++;

  // bottom side
  data.face_colors[tri] = glm::vec4(0.5f, 0.5f, 1.0f, 0);  triangles[tri].x = 0; triangles[tri].y = 4; triangles[tri].z = 1; tri++;
  data.face_colors[tri] = glm::vec4(0.5f, 0.5f, 1.0f, 0);  triangles[tri].x = 1; triangles[tri].y = 4; triangles[tri].z = 5; tri++;

  // top side
  data.face_colors[tri] = glm::vec4(1.0f, 1.0f, 1.0f, 0);  triangles[tri].x = 2; triangles[tri].y = 3; triangles[tri].z = 6; tri++;
  data.face_colors[tri] = glm::vec4(1.0f, 1.0f, 1.0f, 0);  triangles[tri].x = 3; triangles[tri].y = 7; triangles[tri].z = 6; tri++;

  // front side
  data.face_colors[tri] = glm::vec4(0,0,1, 0); triangles[tri].x = 0; triangles[tri].y = 2; triangles[tri].z = 4; tri++;
  data.face_colors[tri] = glm::vec4(0,0,1, 0); triangles[tri].x = 2; triangles[tri].y = 6; triangles[tri].z = 4; tri++;

  // back side
  data.face_colors[tri] = glm::vec4(1,1,0, 0); triangles[tri].x = 1; triangles[tri].y = 5; triangles[tri].z = 3; tri++;
  data.face_colors[tri] = glm::vec4(1,1,0, 0); triangles[tri].x = 3; triangles[tri].y = 5; triangles[tri].z = 7; tri++;

  rtcSetGeometryVertexAttributeCount(mesh,1);
  rtcSetSharedGeometryBuffer(mesh,
                             RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE /* type */,
                             0 /* slot */,
                             RTC_FORMAT_FLOAT3,
                             data.vertex_colors,
                             0 /* byteOffset */,
                             sizeof(glm::vec4) /* byteStride */,
                             8 /* itemCount */);

  rtcCommitGeometry(mesh);
  return mesh;
}


RTCGeometry spheresGeometry(RTCDevice device) {
    RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SPHERE_POINT);

    float radius = 0.1;
    glm::vec4* vertices = (glm::vec4*) rtcSetNewGeometryBuffer(mesh,
                                                               RTC_BUFFER_TYPE_VERTEX,
                                                               0 /* slot */,
                                                               RTC_FORMAT_FLOAT4,
                                                               sizeof(glm::vec4) /* byteStride */,
                                                               8 /* itemCount */);
    vertices[0] = glm::vec4{ 1.0,  1.0,  1.0, radius};
    vertices[1] = glm::vec4{ 1.0, -1.0,  1.0, radius};
    vertices[2] = glm::vec4{-1.0,  1.0,  1.0, radius};
    vertices[3] = glm::vec4{-1.0, -1.0,  1.0, radius};
    vertices[4] = glm::vec4{ 1.0,  1.0, -1.0, radius};
    vertices[5] = glm::vec4{ 1.0, -1.0, -1.0, radius};
    vertices[6] = glm::vec4{-1.0,  1.0, -1.0, radius};
    vertices[7] = glm::vec4{-1.0, -1.0, -1.0, radius};
    rtcCommitGeometry(mesh);

    // Define vertex colors (available later by id hit.primID)
//    data.vertex_colors[0] = glm::vec4(0,0,0, 0);
//    data.vertex_colors[1] = glm::vec4(0,0,1, 0);
//    data.vertex_colors[2] = glm::vec4(0,1,0, 0);
//    data.vertex_colors[3] = glm::vec4(0,1,1, 0);
//    data.vertex_colors[4] = glm::vec4(1,0,0, 0);
//    data.vertex_colors[5] = glm::vec4(1,0,1, 0);
//    data.vertex_colors[6] = glm::vec4(1,1,0, 0);
//    data.vertex_colors[7] = glm::vec4(1,1,1, 0);

    return mesh;
}



std::pair<RTCGeometry, DataVertex> axisSpheresGeometry(RTCDevice device) {
    RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SPHERE_POINT);

    std::size_t N = 8;
    glm::vec4* vertices = (glm::vec4*) rtcSetNewGeometryBuffer(mesh,
                                                               RTC_BUFFER_TYPE_VERTEX,
                                                               0 /* slot */,
                                                               RTC_FORMAT_FLOAT4,
                                                               sizeof(glm::vec4) /* byteStride */,
                                                               N /* itemCount */);
    float radius = 0.1;
    vertices[0] = glm::vec4{ 0.0,  0.0,  -100.0, 100.0}; // ground
    vertices[1] = glm::vec4{ 1.0, 0.0,  0.0, radius};
    vertices[2] = glm::vec4{0.0,  1.0,  0.0, radius};
    vertices[3] = glm::vec4{0.0, 0.0,  1.0, radius};
    vertices[4] = glm::vec4{ 0.0,  0.0, 0.0, radius};
    vertices[5] = glm::vec4{ 0.0,  0.0, 0.0, radius};
    vertices[6] = glm::vec4{ 0.0,  0.0, 0.0, radius};
    vertices[7] = glm::vec4{ 0.0,  0.0, 0.0, radius};
    rtcCommitGeometry(mesh);

    DataVertex data{N};
    // Define vertex colors (available later by id hit.primID)
    data.vertex[0] = glm::vec4(1,1,1, 0);
    data.vertex[1] = glm::vec4(1,0,0, 0);
    data.vertex[2] = glm::vec4(0,1,0, 0);
    data.vertex[3] = glm::vec4(0,0,1, 0);
    data.vertex[4] = glm::vec4(0,0,0, 0);
    data.vertex[5] = glm::vec4(0,0,0, 0);
    data.vertex[6] = glm::vec4(0,0,0, 0);
    data.vertex[7] = glm::vec4(0,0,0, 0);

    return std::make_pair(mesh, std::move(data));
}



std::pair<RTCGeometry, VectorMaterial> spheresSmallGeometry(RTCDevice device) {
    RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SPHERE_POINT);

    std::size_t N = 8;
    glm::vec4* vertices = (glm::vec4*) rtcSetNewGeometryBuffer(mesh,
                                                               RTC_BUFFER_TYPE_VERTEX,
                                                               0 /* slot */,
                                                               RTC_FORMAT_FLOAT4,
                                                               sizeof(glm::vec4) /* byteStride */,
                                                               N /* itemCount */);
    float radius = 1.0;
    vertices[0] = glm::vec4{ 0.0,  0.0,  -100.0, 100.0}; // ground
    vertices[1] = glm::vec4{ -2.0, 0.0,  1.0, radius};
    vertices[2] = glm::vec4{0.0,  0.0,  1.0, radius};
    vertices[3] = glm::vec4{2.0, 0.0,  1.0, radius};
    float radius2 = 0.2;
    vertices[4] = glm::vec4{ -2*radius2,  2.0, radius2, radius2};
    vertices[5] = glm::vec4{ 0.0,  2.0, radius2, radius2};
    vertices[6] = glm::vec4{ 2*radius2,  2.0, radius2, radius2};
    vertices[7] = glm::vec4{ 4*radius2,  2.0, radius2, radius2};
    rtcCommitGeometry(mesh);

    VectorMaterial data{N, Material{MaterialType::Lambertian}};
    data[0].albedo = glm::vec3(0.5,0.5,0.5);

    data[1].type = MaterialType::Dielectric;
    data[1].albedo = glm::vec3(0.7,0.1,0.1); // TODO: not used?
    data[1].ref_idx = 1.5;

    data[2].albedo = glm::vec3(0.1,0.7,0.1);

    data[3].type = MaterialType::Metal;
    data[3].albedo = glm::vec3(0.9,0.9,0.6);
    data[3].fuzz = 0.0f;

    data[4].albedo = glm::vec3(0.7,0.7,0.1);
    data[5].albedo = glm::vec3(0.1,0.7,0.7);
    data[6].albedo = glm::vec3(0.7,0.1,0.7);
    data[7].albedo = glm::vec3(0.1,0.1,0.1);

    return std::make_pair(mesh, std::move(data));
}


RTCScene spheresScene(RTCDevice device) {
    RTCScene cube_scene = rtcNewScene(device);
    if (!cube_scene) { throw "rtcNewScene"; }
    RTCGeometry geometry = spheresGeometry(device);
    /* cube_geometry_id = */ rtcAttachGeometry(cube_scene, geometry);
    rtcReleaseGeometry(geometry);
    rtcCommitScene (cube_scene);
    return cube_scene;
}

std::tuple<RTCScene, unsigned int, VectorMaterial> spheresSmallScene(RTCDevice device) {
    RTCScene scene = rtcNewScene(device);
    if (!scene) { throw "rtcNewScene"; }
    auto [geometry, data] = spheresSmallGeometry(device);
    unsigned int geometry_id = rtcAttachGeometry(scene, geometry);
    rtcReleaseGeometry(geometry);
    rtcCommitScene (scene);
    return std::make_tuple(scene, geometry_id, std::move(data));
}


RTCScene cubeScene(RTCDevice device) {
    RTCScene cube_scene = rtcNewScene(device);
    if (!cube_scene) { throw "rtcNewScene"; }
    RTCGeometry cube_geometry = cubeGeometry(device);
    rtcAttachGeometry(cube_scene, cube_geometry);
    rtcReleaseGeometry(cube_geometry);
    rtcCommitScene (cube_scene);
    return cube_scene;
}


RTCScene cubeX2InstanceScene(RTCDevice device) {
    auto cube_scene = cubeScene(device);

    auto generate_instance = [&](const glm::vec3 translate) -> RTCGeometry{
        RTCGeometry instance0 = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_INSTANCE);
        rtcSetGeometryInstancedScene(instance0, cube_scene); // assign instance0 <- cube_scene
        //glm::mat4 mat = glm::rotate(glm::mat4{1.0}, glm::radians(15.0f), {0.0, 1.0, 0.0});
        glm::mat4 mat = glm::mat4{1.0};
        mat = glm::scale(mat, glm::vec3{0.3f});
        mat = glm::translate(mat, translate);
        rtcSetGeometryTransform(instance0, 0, RTC_FORMAT_FLOAT4X4_COLUMN_MAJOR, (float*)&mat);
        rtcCommitGeometry(instance0);
        return instance0;
    };

    // Main scence with instances
    RTCScene scene = rtcNewScene(device);

    auto instance0 = generate_instance(glm::vec3{0.0, 1.5, 0.0});
    rtcAttachGeometry(scene,instance0);
    rtcReleaseGeometry(instance0);

    auto instance1 = generate_instance(glm::vec3{0.0, -1.5, 0.0});
    rtcAttachGeometry(scene,instance1);
    rtcReleaseGeometry(instance1);

    rtcCommitScene(scene);
    return scene;
}




void testImage() {
    Image im{200,100};
    im.setColor(0,0, glm::vec3{1.0});
    im.setColor(1,1, glm::vec3{1.0, 1.0, 0.0});
    im.setColor(200-1,0, glm::vec3{1.0, 0.0, 0.0});
    im.setColor(0,100-1, glm::vec3{0.0, 1.0, 0.0});
    im.setColor(200-1,100-1, glm::vec3{0.0, 0.0, 1.0});
    im.writeToFileBMP("output.bmp");
}






struct Renderer {
    Image image{frame_columns, frame_rows};
    const float aspect = float{frame_columns} / frame_rows;
    CameraRegular camera{{-3,7,2}, {0,0,0}, {0,0,1}, 40, aspect};

    RTCDevice device;
    RTCScene scene; // TODO: Do I need to keep it here?
    std::unique_ptr<RayTracer> rt;

    glm::vec3 getSample(float xf_rand, const float yf_rand)
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
                image.setColor(xn, frame_rows - 1 - yn, color);
            }
        }
    }

    void renderTile(std::size_t tile) {
        assert(tile < tiles_rows*tiles_columns);
        const std::size_t tile_column = tile % tiles_columns;
        const std::size_t tile_row = tile / tiles_columns;

        const std::size_t c_start = tile_column * TILE_SIZE;
        const std::size_t r_start = tile_row * TILE_SIZE;

        const std::size_t c_stop = c_start + TILE_SIZE;
        const std::size_t r_stop = r_start + TILE_SIZE;
        renderRegion(c_start, c_stop, r_start, r_stop);
    }

    void render() {
        #pragma omp parallel for
        for (std::size_t tile = 0; tile < tiles_rows*tiles_columns; tile++) {
            renderTile(tile);
        }

        // Different experiments with parallelisation
        //    std::vector<std::size_t> tiles_list(tiles_rows*tiles_columns);
        //    for (std::size_t i = 0; i < tiles_list.size(); i++) {
        //        tiles_list[i] = i;
        //    }
        //tbb::parallel_for(std::size_t(0), std::size_t(tiles_rows*tiles_columns), )
        //std::for_each(std::execution::par, std::begin(tiles_list),  std::end(tiles_list), render_tile);

        //image.writeToFileBMP("/home/universe.dart.spb/pwoloszkiewicz/tmp/output.bmp");
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
        this->scene = scene;
        rt = std::make_unique<RayTracer>(scene, std::move(data));
    }

    ~Renderer() {
        // TODO: scene release
        // TODO: free buffers
        rtcReleaseScene(scene);
        rtcReleaseDevice(device);
    }
};

//=========================================================
#include "renderer.h"
namespace DRenderer {

Renderer* create() {
    return new Renderer{};
}

void destroy(Renderer* r) {
    delete r;
}

int numberOfTiles() {
    return tiles_rows*tiles_columns;
}

void renderTile(Renderer* r, int tile_id) {
    r->renderTile(tile_id);
}

void render(Renderer* r) {
    r->render();
}

void writeImage(Renderer* r) {
    r->writeImage("output.bmp");
}

} // namespace: DRenderer
