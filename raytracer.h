#pragma once
#include <embree3/rtcore.h>
#include <vector>
#include "utils.h"


RTCRayHit initRayHit(const glm::vec3& origin = glm::vec3{},
                     const glm::vec3& direction = glm::vec3{});


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


class RayTracer {
    RTCScene scene;
    VectorMaterial data;

public:
    RayTracer(RTCScene scene, VectorMaterial&& data):
        scene{scene},
        data{std::move(data)}
    {
    }

    glm::vec3 rayColor(RTCRayHit& rayh, int depth = 1);
};



