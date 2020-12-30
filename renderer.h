#pragma once

struct Renderer;

class DRenderer {
public:
    DRenderer();
    ~DRenderer();
    static int numberOfTiles();
    void renderAllTiles();
    void renderTile(int tile_id);
    void writeImage();
private:
    Renderer* r;
};
