#pragma once
#include "utils.h"


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
