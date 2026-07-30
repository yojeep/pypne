#include "gtl/phmap.hpp"
#include "voxelmap.h"
#include <iostream>
#include <vector>
// voxel_map.h
#include "gtl/phmap.hpp"
#include <cstddef>
#include <stdexcept>

template <typename KT, typename VT> class VoxelMap {
  // ========== 内部实现，外部完全不可见 ==========
  KT NX_, NY_, NZ_;
  gtl::flat_hash_map<KT, VT> map_;

  size_t encode(KT z, KT y, KT x) const noexcept {
    return static_cast<size_t>(z) * static_cast<size_t>(NY_) * NX_ +
           static_cast<size_t>(y) * static_cast<size_t>(NX_) +
           static_cast<size_t>(x);
  }

public:
  explicit VoxelMap(KT nx, KT ny, KT nz, size_t reserve_n = 0)
      : NX_(nx), NY_(ny), NZ_(nz) {
    if (reserve_n)
      map_.reserve(reserve_n);
  }

  // ---------- 语法 1：函数调用风格（最常用）----------
  // vox(3, 7, 42) = 1.23;
  // double v = vox(3, 7, 42);
  VT &operator()(KT z, KT y, KT x) { return map_[encode(z, y, x)]; }

  // ---------- 语法 2：花括号聚合风格 ----------
  // vox[{3, 7, 42}] = 1.23;
  struct Coord {
    KT z, y, x;
  };
  VT &operator[](Coord c) { return map_[encode(c.z, c.y, c.x)]; }

  // ---------- 语法 3：显式命名接口 ----------
  void emplace(KT z, KT y, KT x, VT val) { map_.emplace(encode(z, y, x), val); }

  [[nodiscard]] VT at(KT z, KT y, KT x) const {
    auto it = map_.find(encode(z, y, x));
    if (it == map_.end())
      throw std::out_of_range("VoxelMap::at: key not found");
    return it->second;
  }

  [[nodiscard]] const VT *find(KT z, KT y, KT x) const noexcept {
    auto it = map_.find(encode(z, y, x));
    return it == map_.end() ? nullptr : &it->second;
  }

  [[nodiscard]] bool contains(KT z, KT y, KT x) const noexcept {
    return map_.contains(encode(z, y, x));
  }

  [[nodiscard]] size_t size() const noexcept { return map_.size(); }
  void clear() noexcept { map_.clear(); }
};
// Key 独立定义
struct Key {
  int z, y, x;

  bool operator==(const Key &o) const {
    return z == o.z && y == o.y && x == o.x;
  }
};

// 哈希器引用 Key
struct GridHasher {
  int NX, NY, NZ;
  GridHasher(int nx, int ny, int nz) : NX(nx), NY(ny), NZ(nz) {}
  size_t operator()(const Key &key) const {
    return key.z * NY * NX + key.y * NX + key.x;
  }
};

int main() {
  auto data = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9};
  auto view = MortonSpan<int, 2>(data.data(), 3, 3);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      std::cout << view[i, j] << std::endl;

  int NX = 100, NY = 200, NZ = 50;
  auto encode = [NX, NY](int z, int y, int x) -> size_t {
    return static_cast<size_t>(z) * NY * NX + static_cast<size_t>(y) * NX +
           static_cast<size_t>(x);
  };

  VoxelMap<int, int> myMap(NX, NY, NZ, 1000000000);
  int count{0};
  int res{0};

  for (int z = 0; z < 100'000; ++z)
    for (int y = 0; y < 100; ++y)
      for (int x = 0; x < 100; ++x)
        if (z != y && y != x && x != z) {
          // myMap.emplace(Key{z, y, x}, 0.0);
          // myMap.emplace(z, y, x, 0);
          res += z * NY * NX + y * NX + x;
          data.emplace_back(0);
          count++;
        }
  std::cout << count << std::endl;
  std::cout << res << std::endl;
  return 0;
}