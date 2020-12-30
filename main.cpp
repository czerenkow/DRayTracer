#include <iostream>
#include "renderer.h"

int main() {
    Renderer *r = DRenderer::create();
    DRenderer::render(r);
    DRenderer::destroy(r);
    return 0;
}

