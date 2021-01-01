#pragma once
#include <vector>
#include <glm/vec3.hpp>


class Image {
public:
    using vec3 = glm::vec3;

    Image(const std::size_t columns, const std::size_t rows):
        columns{columns},
        rows{rows},
        data(columns * rows)
    {
    }

    // Load PFM file
    Image(const char filename[]);

    /// (0,0) -> left top corner
    void setColor(const std::size_t x, const std::size_t y, const vec3& c)
    {
        assert(x < columns);
        assert(y < rows);
        data[x + columns * y] = c;
    }

    vec3 getColor(const std::size_t x, const std::size_t y) const {
        return data[x + columns * y];
    }

    int writeToFileBMP(const char filename[]) const;
    int writeToFilePFM(const char filename[]) const;

    // TODO: create getters (those should not be modified outside)
    std::size_t columns;
    std::size_t rows;
private:
    std::vector<vec3> data;
};
