#pragma once
#include <ostream>
#include <iostream>
#include <embree3/rtcore.h>
#include <vector>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#define GLM_GTX_transform
#include <glm/gtx/transform.hpp>
#define GLM_GTC_constants
#include <glm/gtc/constants.hpp>


inline std::ostream& operator<<(std::ostream& os, const glm::vec3& v) {
    os << '[' << v.x << ',' << v.y << ',' << v.z << ']';
    return os;
}





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

//==================================================
// Materials structs
// TODO: It should be defined somewhere else!
//==================================================


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





