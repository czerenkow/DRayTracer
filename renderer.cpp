#include <iostream>
#include <cstdlib>
#include <embree3/rtcore.h>
#include <embree3/rtcore_common.h>

#include <limits>
#include <tuple>
#include <chrono>
#include <random>
#include <algorithm>
#include <memory>

//#include <execution>
//#include <tbb/tbb.h>

#include "image.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#define GLM_GTX_transform
#include <glm/gtx/transform.hpp>
#define GLM_GTC_constants
#include <glm/gtc/constants.hpp>

using std::cout;
using std::endl;

const float infinity = std::numeric_limits<float>::infinity();

const std::size_t TILE_SIZE = 128;
const std::size_t frame_columns_exp = 2600;
const std::size_t frame_rows_exp = 2000;
const int depth_max = 50;
const int samples = 500;

// Calculate colums x rows such that were divided by TILE_SIZE
const std::size_t frame_columns = frame_columns_exp / TILE_SIZE * TILE_SIZE;
const std::size_t frame_rows = frame_rows_exp / TILE_SIZE * TILE_SIZE;
const std::size_t tiles_columns = frame_columns / TILE_SIZE;
const std::size_t tiles_rows = frame_rows / TILE_SIZE;


std::ostream& operator<<(std::ostream& os, const glm::vec3& v) {
    os << '[' << v.x << ',' << v.y << ',' << v.z << ']';
    return os;
}


class RandomGenerator {
public:
    RandomGenerator() {
        //displaySeed();
    }

//    RandomGenerator(int seed): seed{seed} {
//        //displaySeed();
//    }

    void displaySeed() {
        std::cout << "Seed: " << this->seed << '\n';
    }

    float get1() {
        return distribution1(generator);
    }
    float get01() {
        return distribution01(generator);
    }
    int get_int(int i) {
        return dist_int(generator) % i;
    }
private:
    const int seed = std::chrono::system_clock::now().time_since_epoch().count();
    //const long seed = 0; // Seed for testing
    std::default_random_engine generator;
    //static std::uniform_real_distribution<double> distribution05(-0.5f, 0.5f);
    std::uniform_real_distribution<float> distribution1{-1.0, 1.0};
    std::uniform_real_distribution<float> distribution01{0.0, 1.0};
    std::uniform_int_distribution<int> dist_int{0, 1000000}; // TODO: find something better
};

thread_local RandomGenerator rnd;


glm::vec3 random_cos_direction() {
    const float r1 = rnd.get01();
    const float r2 = rnd.get01();
    const float phi = 2.0 * r1 * glm::pi<float>();
    const float sqrt_2r2 = std::sqrt(r2); // In the book, here there is *2 additionally, and as I understand this is incorrect.
    const float x = std::cos(phi) * sqrt_2r2;
    const float y = std::sin(phi) * sqrt_2r2;
    const float z = std::sqrt(1.0 - r2);
    return glm::vec3{x, y, z};
}

/// n must be normal
/// Ortogonal matrix (0,0,1) -> n
///  TEST: m * glm::transpose(m) == I
glm::mat3 build_onb_matrix(const glm::vec3& n) {
    const glm::vec3& axis2 = n;
    const glm::vec3 a = std::abs(axis2.x) > 0.9f ? glm::vec3{0.0, 1.0, 0.0} : glm::vec3{1.0, 0.0, 0.0};
    const glm::vec3 axis1 = glm::normalize( glm::cross(axis2, a) );
    const glm::vec3 axis0 = glm::cross(axis2, axis1);
    return glm::mat3{axis0, axis1, axis2};
}

// normal: must be normalized
// scattered_direction: must be normalized
float lambertian_scattering_pdf(const glm::vec3& normal, const glm::vec3& scattered_direction) {
    const float cosine = glm::dot(normal, scattered_direction);
    return cosine < 0.0 ? 0.0 : cosine / glm::pi<double>();
}



inline glm::vec3 reflect(const glm::vec3& v, const glm::vec3& n) {
    return v - 2.0f * glm::dot(v,n)*n;
}


inline glm::vec3 random_in_unit_sphere() {
    for(;;) {
        glm::vec3 result {rnd.get1(), rnd.get1(), rnd.get1()};
        if (glm::dot(result, result) < 1.0) {
            return result;
        }
    }
}

// ray_in_direction: must be normalized
// normal: must be normalized
bool metal_scatter(const glm::vec3& ray_in_direction, // incoming ray (normalized)
                   const glm::vec3& normal, // surfice normal (normalized)
                   float fuzz,
                   glm::vec3& scattered_dir) {
    /*
     * This:
     *   glm::dot(ray_in_normal, normal) > 0.0
     * means that the ray hit wrong side of the surfice.
     * Even if we check it here it will not guarantee that with 'fuzz'
     * we will get direction outside the surfice, so we need to check it
     * on the end of computing.
     */
    const glm::vec3 reflected = reflect(ray_in_direction, normal);
    scattered_dir = reflected + fuzz * random_in_unit_sphere();
    return glm::dot(scattered_dir, normal) > 0.0;
}



static bool refract(const glm::vec3& v, const glm::vec3& n, float ni_over_nt, glm::vec3& refracted) {
    const glm::vec3 uv = glm::normalize(v);
    const float dt = glm::dot(uv, n);
    const float discriminant = 1.0 - ni_over_nt*ni_over_nt*(1.0 - dt*dt);
    if (discriminant > 0.0) {
        refracted = ni_over_nt * (uv - n*dt) - n*std::sqrt(discriminant);
        return true;
    }
    return false;
}

static double schlick(float cosine, float ref_idx) {
    float r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * std::pow(1.0 - cosine, 5.0);
}



bool dielectric_scatter(const glm::vec3& ray_in_direction, // incoming ray
                        const glm::vec3& normal, // surfice normal
                        float ref_idx,
                        glm::vec3& scattered)
{
    glm::vec3 outward_normal;
    float ni_over_nt;
    float cosine;
    if (glm::dot(ray_in_direction, normal) > 0) {
        outward_normal = -normal;
        ni_over_nt = ref_idx;
        cosine = ref_idx * dot(ray_in_direction, normal) / glm::length(ray_in_direction);
    } else {
        outward_normal = normal;
        ni_over_nt = 1.0 / ref_idx;
        cosine =          -dot(ray_in_direction, normal) / glm::length(ray_in_direction);
    }

    glm::vec3 refracted;
    float reflect_prob;
    if ( refract(ray_in_direction, outward_normal, ni_over_nt, refracted) ) {
        reflect_prob = schlick(cosine, ref_idx);
    } else {
        reflect_prob = 1.0;
    }

    if (rnd.get01() < reflect_prob) {
        scattered = reflect(ray_in_direction, normal);
    } else {
        scattered = refracted;
    }
    return true;
}


class Pdf {
public:
    /**
     * @brief generate generates random vector
     * @return
     */
    virtual glm::vec3 generate() const = 0;

    /**
     * @brief value probability of vector 'direction'
     * @param direction must be normalized
     * @return
     */
    virtual float value(const glm::vec3& direction) const = 0;

    virtual ~Pdf() {}
};


class PdfConsine: public Pdf {
public:
    /**
     * @brief pdf_consine generates cosine distribution of directions on the hemisphere in direction w
     * @param w
     */
    PdfConsine(const glm::vec3& w):
        onb {build_onb_matrix(w)}
    {
    }

    // direction: must be normalized
    float value(const glm::vec3& direction) const override {
        const glm::vec3 onb_main_direction = onb[2]; // onb * [0,0,1]^T    TODO: copy??
        const float cosine = glm::dot(direction, onb_main_direction);
        return cosine <= 0.0f ? 0.0f : cosine / glm::pi<float>();
    }

    glm::vec3 generate() const override {
        // TODO: Is it already normal?
        return onb * random_cos_direction();
    }
private:
    const glm::mat3 onb;
};



enum struct MaterialType {
    Lambertian, Metal, Dielectric, DiffuseLight
};

struct Material {
    MaterialType type {};

    glm::vec3 emmited{}; // DiffuseLight
    glm::vec3 albedo{}; // all
    float fuzz{}; // Metal
    float ref_idx{}; // Dielectric

    Material(MaterialType type): type{type} {}
    Material() {}
};


using VectorMaterial = std::vector<Material>;


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



inline glm::vec3 rayDirection(const RTCRay& rayh) {
    return {rayh.dir_x, rayh.dir_y, rayh.dir_z};
}

inline glm::vec3 rayOrigin(const RTCRay& rayh) {
    return {rayh.org_x, rayh.org_y, rayh.org_z};
}

inline void raySetOrigin(RTCRay& ray, const glm::vec3& origin) {
    ray.org_x = origin.x;
    ray.org_y = origin.y;
    ray.org_z = origin.z;
}

inline void raySetDirection(RTCRay& ray, const glm::vec3& dir) {
    ray.dir_x = dir.x;
    ray.dir_y = dir.y;
    ray.dir_z = dir.z;
}


inline glm::vec3 hitNormal(const RTCHit& rayh) {
    return {rayh.Ng_x, rayh.Ng_y, rayh.Ng_z};
}

//inline glm::vec3 hitPoint(const RTCRay& rayh) {
//    return rayOrigin(rayh) + rayh.tfar * rayDirection(rayh);
//}


RTCRayHit initRayHit(const glm::vec3& origin = glm::vec3{},
                     const glm::vec3& direction = glm::vec3{})
{
    RTCRayHit rhit;
    RTCRay& ray = rhit.ray;
    ray.tnear = 0.0001f; // acne problem
    ray.tfar = infinity;

    ray.flags = 0;
    ray.mask = 0;

    ray.time = 0.0f;
    ray.id = 0;

    rhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    raySetOrigin(rhit.ray, origin);
    raySetDirection(rhit.ray, direction);
    return rhit;
}


/**
 * @brief getSkyColor color of space around the space.
 * @param r
 * @return
 */
glm::vec3 getBackgroundColor(const RTCRay& r) {
    glm::vec3 unit_direction = glm::normalize(rayDirection(r));
    float t = 0.5*(unit_direction.z + 1.0);
    return (1.0f - t)*glm::vec3(1.0) + t*glm::vec3(0.5, 0.7, 1.0);
}


class RayTracer {
    RTCScene scene;
    VectorMaterial data;

public:
    RayTracer(RTCScene scene, VectorMaterial&& data):
        scene{scene},
        data{std::move(data)}
    {
    }

    glm::vec3 rayColor(RTCRayHit& rayh, int depth = 1) {
        if (depth == depth_max) {
            return glm::vec3{0.0, 0.0, 0.0};
        }

        RTCIntersectContext context;
        rtcInitIntersectContext(&context);
        rtcIntersect1(scene, &context, &rayh);
        RTCHit& hit = rayh.hit;
        if (hit.geomID == RTC_INVALID_GEOMETRY_ID) {
            return getBackgroundColor(rayh.ray);
        }
        Material& material = data[hit.primID];

        RTCGeometry geometry = rtcGetGeometry(scene, hit.geomID);
        glm::vec4* vertex_data = (glm::vec4*)rtcGetGeometryBufferData(geometry, RTC_BUFFER_TYPE_VERTEX, 0 /*slot*/);
        glm::vec4& vertex_obj = vertex_data[hit.primID];
        const glm::vec3 hit_normal_n = hitNormal(rayh.hit);
        const glm::vec3 hit_normal = glm::normalize(hit_normal_n);
        const glm::vec3 hit_point = glm::vec3{vertex_obj.x, vertex_obj.y, vertex_obj.z} + hit_normal_n;

        switch (material.type) {
        case MaterialType::Lambertian: {
            const PdfConsine p0{hit_normal}; // Cosine (lambertian) sampling
            const glm::vec3 scattered_ray_dir = glm::normalize(p0.generate()); // TODO: do we need to normalize here?
            const float pdf = p0.value(scattered_ray_dir);
            const float scattering_pdf = lambertian_scattering_pdf(hit_normal, scattered_ray_dir);

            RTCRayHit scattered_ray = initRayHit(hit_point, scattered_ray_dir);
            auto c = rayColor(scattered_ray, depth + 1);
            return material.albedo * c * scattering_pdf / pdf;
        }
        case MaterialType::Metal: {
            const glm::vec3 ray_dir = glm::normalize(rayDirection(rayh.ray));
            glm::vec3 scattered_ray_dir;
            // TODO: Calculation with fuzz>0 is incorrect. This is visible on the edge of the ball.
            if (!metal_scatter(ray_dir, hit_normal, material.fuzz, scattered_ray_dir)) {
                // ERROR
                return glm::vec3{1,0,0}; // TODO: why emmited?
            }

            // Issue next ray
            RTCRayHit scattered_ray = initRayHit(hit_point, scattered_ray_dir);
            auto c = rayColor(scattered_ray, depth + 1);

            return material.albedo * c;
        }
        case MaterialType::Dielectric: {
            const glm::vec3 ray_dir = glm::normalize(rayDirection(rayh.ray));
            glm::vec3 scattered_ray_dir;
            if (!dielectric_scatter(ray_dir, hit_normal, material.ref_idx, scattered_ray_dir)) {
                // ERROR
                return material.emmited; // TODO: why emmited?
            }

            // Issue next ray
            RTCRayHit scattered_ray = initRayHit(hit_point, scattered_ray_dir);
            auto c = rayColor(scattered_ray, depth + 1);
            return /*material.albedo **/ c;
        }
        default:
            throw "Unknown material";
        }
    }

};






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
                    const float yf_rand = yf + rnd.get01();
                    const float xf_rand = xf + rnd.get01();
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
