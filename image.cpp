#include "image.h"
#include <fstream>
#include <ios>
#include <iostream>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"


int Image::writeToFileBMP(const char filename[]) const
{
    // convert to 8 bit/channel color
    struct vec3_uint8 {
        uint8_t r,g,b;
    };
    std::vector<vec3_uint8> data_out(data.size());
    std::transform(data.cbegin(), data.cend(), data_out.begin(), [](const vec3& c) -> vec3_uint8 {
        using VT = vec3::value_type;
        assert(VT{0} <= c.r && c.r <= VT{1});
        assert(VT{0} <= c.g && c.g <= VT{1});
        assert(VT{0} <= c.b && c.b <= VT{1});
        auto c2 = VT{255} * c;
        return {uint8_t(c2.r), uint8_t(c2.g), uint8_t(c2.b)};
    });

    return stbi_write_bmp(filename, static_cast<int>(columns), static_cast<int>(rows), 3, data_out.data());
}


int Image::writeToFilePFM(const char filename[]) const {
    std::ofstream file{filename, std::ios::out | std::ios::binary};
    if (!file.is_open()) {
        std::cerr << "ERROR: writeToFilePFM: open file\n";
        return 0;
    }

    file << "PF" << '\n';
    file << columns << ' ' << rows << '\n';
    file << "-1.0" << '\n';

    const std::size_t row_size = columns * sizeof(vec3);
    for (std::size_t y = rows; y >= 1; y--) {
        std::size_t y_curr = y - 1;
        // little endian
        auto row_begin = reinterpret_cast<const char*>(&data[y_curr * columns]);
        file.write(row_begin, static_cast<std::streamsize>(row_size));
    }
    return 1;
}


// Load PFM file
Image::Image(const char filename[]) {
    std::fstream file(filename);
    if (!file.is_open()) {
        throw "image::image(): can not open file";
    }

    std::string line;
    file >> line;
    if (line != "PF") {
        throw "not PMF file?";
    }
    file >> columns;
    file >> rows;

    std::cout << "Loading PFM file. Size: " << columns << 'x' << rows << '\n';
    file >> line;
    if (line != "-1.0") {
        throw "PMF file unsupported endianess";
    }
    const std::size_t elem_size = sizeof(vec3);
    const std::size_t row_size = columns * elem_size;

    file.get();// read '\n'
    data.resize(rows * columns);
    for (std::size_t y = rows; y >= 1; y--) {
        std::size_t y_curr = y - 1;
        // little endian
        file.read(reinterpret_cast<char*>(&data[y_curr * columns]), row_size);
    }
}

