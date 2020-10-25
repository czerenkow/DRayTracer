#include <iostream>
#include "renderer.h"

int main() {
//    cout << glm::dot(glm::vec3{0.13, -0.98, 0.073}, glm::vec3{0.11, 0.36, -0.925}) << endl;
//    cout << glm::dot(-glm::vec3{0.219, -0.69, 0.69}, glm::vec3{0.11, 0.36, -0.925}) << endl;
    Renderer *r = DRenderer::create();
    DRenderer::render(r);
    DRenderer::destroy(r);
    return 0;
}

