#pragma once

class Renderer;

namespace DRenderer {

Renderer* create();
void destroy(Renderer* r);
int numberOfTiles();
void render(Renderer* r);
void renderTile(Renderer* r, int tile_id);
void writeImage(Renderer* r);

}
