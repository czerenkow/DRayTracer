#pragma once
#include "constants.h"
#include <cassert>

struct TileInfo {
    TileInfo(std::size_t tile) {
        assert(tile < numberOfTiles());
        tile_column = tile % tiles_columns;
        tile_row = tile / tiles_columns;

        c_start = tile_column * TILE_SIZE;
        r_start = tile_row * TILE_SIZE;

        c_stop = c_start + TILE_SIZE;
        r_stop = r_start + TILE_SIZE;
    }

    static std::size_t numberOfTiles() {
        return tiles_rows * tiles_columns;
    }

    std::size_t tile_column;
    std::size_t tile_row;

    std::size_t c_start;
    std::size_t r_start;

    std::size_t c_stop;
    std::size_t r_stop;
};



struct Renderer;
class Image;

class DRenderer {
public:
    DRenderer();
    ~DRenderer();
    void renderAllTiles();
    void renderTile(int tile_id);
    Image& getImageRef();
    void writeImage();
private:
    Renderer* r;
};


void testImage();
