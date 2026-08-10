#pragma once

#include <memory>
#include <thread_pool/BS_thread_pool.hpp>
#include <Eigen/CXX11/Tensor>
inline int num_workers(1);
class GlobalThreadPool {
private:
  // 使用 inline 让静态变量可在头文件中定义
  inline static std::unique_ptr<BS::thread_pool<>> instance = nullptr;
  inline static std::once_flag init_flag;

public:
  static BS::thread_pool<> &get() {
    // 线程安全的懒加载
    if (!instance) {
      instance = std::make_unique<BS::thread_pool<>>(num_workers);
    }
    return *instance;
  }

  static void shutdown() {
    if (instance) {
      instance.reset(); // 释放资源
    }
  }

  // 禁止实例化
  GlobalThreadPool() = delete;
  GlobalThreadPool(const GlobalThreadPool &) = delete;
  GlobalThreadPool &operator=(const GlobalThreadPool &) = delete;
};

template <typename T> inline T sq(T x) noexcept { return x * x; }
template <typename T> inline T norm_sq(T x, T y, T z) noexcept {
  return x * x + y * y + z * z;
}
template <typename T> inline T norm(T x, T y, T z) noexcept {
  return std::sqrt(norm_sq(x, y, z));
}

template <typename Func>
inline void for_each_voxel_in_sphere_delta(int cz, int cy, int cx, float radius,
                                           int nz, int ny, int nx,
                                           Func &&func) noexcept {
  float r2 = radius * radius;
  int z_start = std::max(static_cast<int>(std::ceil(cz - radius)), 0);
  int z_end = std::min(static_cast<int>(std::floor(cz + radius)) + 1, nz);

  for (int z = z_start; z < z_end; ++z) {
    float dz = z - cz;
    float ry2 = r2 - sq(dz);
    if (ry2 <= 0)
      continue;
    float ry = std::sqrtf(ry2);

    int y_start = std::max(static_cast<int>(std::ceil(cy - ry)), 0);
    int y_end = std::min(static_cast<int>(std::floor(cy + ry)) + 1, ny);

    for (int y = y_start; y < y_end; ++y) {
      float dy = y - cy;
      float rx2 = ry2 - sq(dy);
      if (rx2 <= 0)
        continue;
      float rx = std::sqrtf(rx2);

      int x_start = std::max(static_cast<int>(std::ceil(cx - rx)), 0);
      int x_end = std::min(static_cast<int>(std::floor(cx + rx)) + 1, nx);

      for (int x = x_start; x < x_end; ++x) {
        float dx = x - cx;
        func(z, y, x, dz, dy, dx);
      }
    }
  }
}

template <typename Func>
inline void for_each_voxel_in_sphere(int cz, int cy, int cx, float radius,
                                     int nz, int ny, int nx,
                                     Func &&func) noexcept {
  for_each_voxel_in_sphere_delta(
      cz, cy, cx, radius, nz, ny, nx,
      [&func](int z, int y, int x, float, float, float) { func(z, y, x); });
}
