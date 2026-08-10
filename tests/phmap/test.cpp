#include <Eigen/CXX11/Tensor>
#include <cstddef>
#include <gtl/phmap.hpp>
#include "thread_pool/BS_thread_pool.hpp"
#include <iostream>

inline int num_workers(64);
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

struct vxl {
  vxl() = default;
  vxl(int z_in, int y_in, int x_in) : z(z_in), y(y_in), x(x_in) {}
  // 1. 禁止拷贝
  vxl(const vxl &) = delete;
  vxl &operator=(const vxl &) = delete;

  // 2. 允许移动
  vxl(vxl &&) noexcept = default;
  vxl &operator=(vxl &&) noexcept = default;
  int z;
  int y;
  int x;
  bool operator==(const vxl &other) const noexcept {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct IdentityHash {
  inline std::size_t operator()(int64_t key) const noexcept {
    return static_cast<std::size_t>(key);
  }
};

class VxlHashMap {
public:
  using MapType = gtl::parallel_flat_hash_map<
      int64_t, vxl, IdentityHash, std::equal_to<int64_t>,
      std::allocator<std::pair<const int64_t, vxl>>, 12, std::mutex>;
  MapType map_;
  int nz, ny, nx;

  VxlHashMap(int nz_in, int ny_in, int nx_in)
      : nz(nz_in), ny(ny_in), nx(nx_in) {}

  // 1. 核心转换公式
  inline int64_t to_key(const vxl &v) const noexcept {
    return static_cast<int64_t>(v.z) * nx * ny +
           static_cast<int64_t>(v.y) * nx + v.x;
  }

  // 2. 极简并行插入：只需传入 total，内部自动算坐标并插入
  void parallel_insert(size_t total) {
    auto &pool = GlobalThreadPool::get();
    pool.detach_blocks(size_t(0), total, [&](size_t start, size_t end) {
      for (size_t i = start; i < end; ++i) {
        int x = static_cast<int>(i % nx);
        int y = static_cast<int>((i / nx) % ny);
        int z = static_cast<int>(i / (static_cast<size_t>(nx) * ny));
        vxl v{z, y, x};
        int64_t key = to_key(v);
        map_.lazy_emplace(key,
                          [&](const auto &ctor) { ctor(key, std::move(v)); });
      }
    });
    pool.wait(); // 阻塞等待完成
  }
};

// ================= 测试 =================
int main() {
  int nx = 1000, ny = 1000, nz = 10000;
  VxlHashMap map(nz, ny, nx);

  size_t total = static_cast<size_t>(nx) * ny * nz;
  map.map_.reserve(total);

  // 极其清爽！只传一个 total 即可
  map.parallel_insert(total);

  std::cout << "Map size: " << map.map_.size() << std::endl;
  GlobalThreadPool::shutdown();
  return 0;
};