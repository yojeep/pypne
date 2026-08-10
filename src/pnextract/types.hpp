#pragma once
#include <limits>
#include <numbers>
#include <cstdint>
#include <Eigen/CXX11/Tensor>

namespace Constants {
typedef long long int64_t;
typedef float float32_t;
typedef double float64_t;

constexpr float32_t F32_EPS = std::numeric_limits<float>::epsilon();
constexpr float64_t F64_EPS = std::numeric_limits<double>::epsilon();
constexpr float32_t F32_PI = std::numbers::pi;
constexpr float64_t F64_PI = std::numbers::pi;
constexpr int32_t I32_MAX = std::numeric_limits<int32_t>::max();
constexpr int32_t I32_MIN = std::numeric_limits<int32_t>::min();
constexpr int64_t I64_MAX = std::numeric_limits<int64_t>::max();
constexpr int64_t I64_MIN = std::numeric_limits<int64_t>::min();

static constexpr int NEIGHBOR_OFFSETS_3[3][3] = {
    {0, 0, 1}, // axis=0
    {0, 1, 0}, // axis=1
    {1, 0, 0}  // axis=2
};
} // namespace Constants

typedef Eigen::Array<float, 1, 3, Eigen::RowMajor> Array3f32;
typedef Eigen::Array<float, 1, 2, Eigen::RowMajor> Array2f32;
typedef Eigen::Matrix<int, 1, 3, Eigen::RowMajor> Vector3i32;
typedef Eigen::Matrix<int, 1, 2, Eigen::RowMajor> Vector2i32;
typedef Eigen::Matrix<int64_t, 1, 3, Eigen::RowMajor> Vector3i64;
typedef Eigen::Matrix<int64_t, 1, 2, Eigen::RowMajor> Vector2i64;
typedef Eigen::Matrix<float, 1, 3, Eigen::RowMajor> Vector3f32;
typedef Eigen::Matrix<float, 1, 2, Eigen::RowMajor> Vector2f32;
typedef Eigen::Matrix<double, 1, 3, Eigen::RowMajor> Vector3f64;
typedef Eigen::Matrix<double, 1, 2, Eigen::RowMajor> Vector2f64;
typedef Eigen::Tensor<int, 3, Eigen::RowMajor> TensorXXXDi32;
typedef Eigen::Tensor<float, 3, Eigen::RowMajor> TensorXXXDf32;
