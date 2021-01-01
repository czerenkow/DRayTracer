#pragma once
#include <cstddef>
#include <limits>

inline constexpr int depth_max = 50;

inline constexpr std::size_t TILE_SIZE = 128;
inline constexpr std::size_t frame_columns_exp = 3000;
inline constexpr std::size_t frame_rows_exp = 2000;
inline constexpr int samples = 700;

// Calculate colums x rows such that were divided by TILE_SIZE
inline constexpr std::size_t frame_columns = frame_columns_exp / TILE_SIZE * TILE_SIZE;
inline constexpr std::size_t frame_rows = frame_rows_exp / TILE_SIZE * TILE_SIZE;
inline constexpr std::size_t tiles_columns = frame_columns / TILE_SIZE;
inline constexpr std::size_t tiles_rows = frame_rows / TILE_SIZE;

constexpr float infinity = std::numeric_limits<float>::infinity();
