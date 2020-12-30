#include "scene.h"

//============================================================
// Definition of redenred scene.
// Most of function defined here are not used. Those are
// for testing purposes or experimentations.
//============================================================



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




