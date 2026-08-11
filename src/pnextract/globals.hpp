#pragma once

#include <memory>
#include <thread_pool/BS_thread_pool.hpp>
#include <Eigen/CXX11/Tensor>

inline int num_workers(1);
inline constexpr float _0p5(0.5);

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

template <typename Func, typename Tc, typename Tr>
inline void for_each_voxel_in_sphere(Tc zc, Tc yc, Tc xc, Tr radius, int nz,
                                     int ny, int nx, Func &&func) noexcept {
  Tr r2 = radius * radius;
  int z_start = std::max(static_cast<int>(std::ceil(zc - radius)), 0);
  int z_end = std::min(static_cast<int>(std::floor(zc + radius)) + 1, nz);

  for (int z = z_start; z < z_end; ++z) {
    Tr dz = z - zc;
    Tr ry2 = r2 - sq(dz);
    if (ry2 <= 0)
      continue;
    Tr yr = std::sqrt(ry2);

    int y_start = std::max(static_cast<int>(std::ceil(yc - yr)), 0);
    int y_end = std::min(static_cast<int>(std::floor(yc + yr)) + 1, ny);

    for (int y = y_start; y < y_end; ++y) {
      Tr dy = y - yc;
      Tr rx2 = ry2 - sq(dy);
      if (rx2 <= 0)
        continue;
      Tr xr = std::sqrt(rx2);

      int x_start = std::max(static_cast<int>(std::ceil(xc - xr)), 0);
      int x_end = std::min(static_cast<int>(std::floor(xc + xr)) + 1, nx);

      for (int x = x_start; x < x_end; ++x) {
        func(z, y, x);
      }
    }
  }
}