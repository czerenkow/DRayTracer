#include "raytracer.h"
#include "random.h"
#include "constants.h"
#include "utils.h"



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


RTCRayHit initRayHit(const glm::vec3& origin,
                     const glm::vec3& direction)
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





glm::vec3 random_cos_direction() {
    const float r1 = randomns::get01();
    const float r2 = randomns::get01();
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
        glm::vec3 result {randomns::get1(), randomns::get1(), randomns::get1()};
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

    if (randomns::get01() < reflect_prob) {
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



glm::vec3 RayTracer::rayColor(RTCRayHit& rayh, int depth) {
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
