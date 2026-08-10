#include "../../src/include/bit_vector.hpp"
#include <cstddef>
struct point {
  size_t x, y, z;
};
int main() {
  // 1. 创建 3D 位向量 10x10x10
  pasta::BitVector bv_(10 * 10 * 10, false);
  BitVector3D bv(bv_, 10, 10, 10);
  // 2. 填充数据
  for (size_t x = 0; x < 10; ++x) {
    bv(x, x, x) = true;
  }
  for (size_t ind = 0; ind < bv.size(); ++ind)
    bv[ind] = true;

  // 3. 数据准备完毕，手动调用 build
  bv.buildRankQuery();

  // 4. 直接通过对象调用 rank 查询
  std::cout << "Total 1s in vector: " << bv.rank1(9, 9, 9) + 1
            << std::endl; // 输出 10
  auto const val = bv(0);
  if (bv(0))
    std::cout << "true" << std::endl;
  else
    std::cout << "false" << std::endl;

  std::cout << "Total 0s in vector: " << bv.count0() << std::endl;

  std::cout << bool(bv(1, 1, 1)) << std::endl;

  BitMap3D<point> bm(bv_, 10, 10, 10);
  bm.bv_[0] = true;
  bm.buildRankQuery();
  // bm.vec_data_.reserve(bm.count1());
  bm.vec_data_.emplace_back(99, 99, 99);
  bm.vec_data_.emplace_back(99, 99, 99);

  bm.bv_[1] = true;
  bm.buildRankQuery();
  point &p = bm(0, 0, 1);
  p.y = 100;
  std::cout << bm(0, 0, 1).y << std::endl;

  return 0;
}