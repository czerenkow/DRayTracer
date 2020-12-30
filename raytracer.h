#pragma once
#include <vector>
#include "utils.h"


RTCRayHit initRayHit(const glm::vec3& origin = glm::vec3{},
                     const glm::vec3& direction = glm::vec3{});



class RayTracer {
    RTCScene scene;
    VectorMaterial data;

public:
    RayTracer(RTCScene scene, VectorMaterial&& data):
        scene{scene},
        data{std::move(data)}
    {}

    glm::vec3 rayColor(RTCRayHit& rayh, int depth = 1);
};



