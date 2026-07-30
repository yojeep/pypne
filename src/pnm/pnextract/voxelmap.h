#pragma once
#include "gtl/phmap.hpp"
#include "gtl/phmap_utils.hpp" // minimal header providing phmap::HashState()
#include "libmorton/morton.h"
#include <cstdint>
#include <iostream>
#include <cstddef>
#include <mdspan>

// namespace utils {

// 3. 自定义布局策略：MortonLayout
// struct MortonLayout {
//   template <class Extents> struct mapping {
//     using index_type = typename Extents::index_type; // ← 跟着 Extents 走

//     Extents exts_;

//     constexpr mapping(Extents &e) noexcept : exts_(e) {}
//     constexpr Extents extents() const noexcept { return exts_; }

//     constexpr index_type operator()(index_type y, index_type x) const noexcept {
//       x = static_cast<uint_fast32_t>(x);
//       y = static_cast<uint_fast32_t>(y);
//       auto answer =
//           static_cast<index_type>(libmorton::morton2D_64_encode(y, x));
//           std::cout << answer << std::endl;
//       return answer;
//     }

//     constexpr index_type operator()(index_type z, index_type y,
//                                     index_type x) const noexcept {
//       x = static_cast<uint_fast32_t>(x);
//       y = static_cast<uint_fast32_t>(y);
//       z = static_cast<uint_fast32_t>(z);
//       auto answer =
//           static_cast<index_type>(libmorton::morton3D_64_encode(x, y, z));
//       return answer;
//     }
//   };
// };
// } // namespace utils

// 4. 提供便捷的类型别名（Type Alias）
// 将冗长的模板参数封装，业务代码只需使用 MortonSpan3D 即可
template <typename T, size_t Rank>
using MortonSpan =
    std::mdspan<T, std::dextents<size_t, Rank>>;

// std::vector<int> data = {0, 1, 2, 3, 4, 5, 6, 7, 8};
// auto span = MortonSpan<int, 2>(data.data(), 3, 3);
