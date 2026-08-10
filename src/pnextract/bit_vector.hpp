#pragma once
#include <cstddef>
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/wide_rank.hpp>
#include <vector>

class BitVector3D {
public:
  pasta::BitVector &bv_;
  size_t nz, ny, nx;
  pasta::WideRank<> rank_query_;

public:
  explicit BitVector3D(pasta::BitVector &bv, size_t nz, size_t ny, size_t nx)
      : bv_(bv), nz(nz), ny(ny), nx(nx) {}

  void buildRankQuery() { rank_query_ = pasta::WideRank<>(bv_); }

  size_t size() noexcept { return bv_.size(); }
  size_t size() const noexcept { return bv_.size(); }

  size_t ravel_idx(size_t z, size_t y, size_t x) noexcept {
    return (z * ny + y) * nx + x;
  }
  size_t ravel_idx(size_t z, size_t y, size_t x) const noexcept {
    return (z * ny + y) * nx + x;
  }

  bool isinside(size_t idx) noexcept { return 0 <= idx && idx < size(); }
  bool isinside(size_t idx) const noexcept { return 0 <= idx && idx < size(); }

  bool isinside(size_t z, size_t y, size_t x) noexcept {
    return 0 <= z && z < nz && 0 <= y && y < ny && 0 <= x && x < nx;
  }
  bool isinside(size_t z, size_t y, size_t x) const noexcept {
    return 0 <= z && z < nz && 0 <= y && y < ny && 0 <= x && x < nx;
  }

  pasta::BitAccess operator[](size_t idx) noexcept { return bv_[idx]; }
  bool operator[](size_t idx) const noexcept { return bv_[idx]; }

  pasta::BitAccess operator[](size_t z, size_t y, size_t x) noexcept {
    return bv_[ravel_idx(z, y, x)];
  }
  bool operator[](size_t z, size_t y, size_t x) const noexcept {
    return bv_[ravel_idx(z, y, x)];
  }

  pasta::BitAccess operator()(size_t idx) noexcept { return bv_[idx]; }
  bool operator()(size_t idx) const noexcept { return bv_[idx]; }

  pasta::BitAccess operator()(size_t z, size_t y, size_t x) noexcept {
    return bv_[ravel_idx(z, y, x)];
  }
  bool operator()(size_t z, size_t y, size_t x) const noexcept {
    return bv_[ravel_idx(z, y, x)];
  }

  size_t rank0(size_t idx) noexcept { return rank_query_.rank0(idx); }
  size_t rank0(size_t idx) const noexcept { return rank_query_.rank0(idx); }

  size_t rank0(size_t z, size_t y, size_t x) noexcept {
    return rank_query_.rank0(ravel_idx(z, y, x));
  }
  size_t rank0(size_t z, size_t y, size_t x) const noexcept {
    return rank_query_.rank0(ravel_idx(z, y, x));
  }

  size_t rank1(size_t idx) noexcept { return rank_query_.rank1(idx); }
  size_t rank1(size_t idx) const noexcept { return rank_query_.rank1(idx); }

  size_t rank1(size_t z, size_t y, size_t x) noexcept {
    return rank_query_.rank1(ravel_idx(z, y, x));
  }
  size_t rank1(size_t z, size_t y, size_t x) const noexcept {
    return rank_query_.rank1(ravel_idx(z, y, x));
  }

  size_t count0() noexcept { return rank0(size()); }
  size_t count0() const noexcept { return rank0(size()); }

  size_t count1() noexcept { return rank1(size()); }
  size_t count1() const noexcept { return rank1(size()); }
};

template <typename T> class BitMap3D : public BitVector3D {
public:
  std::vector<T> vec_data_;
  T None;

public:
  explicit BitMap3D(pasta::BitVector &bv, size_t nz, size_t ny, size_t nx)
      : BitVector3D(bv, nz, ny, nx), None(T{}) {}

  // 非 const 版本 → 返回 T&（不变）
  T &operator()(size_t idx, bool check = false) noexcept {
    if (check && !isinside(idx))
      return None;
    if (bv_[idx])
      return vec_data_[rank1(idx)];
    else
      return None;
  }

  // const 版本 → 返回 const T& ← 改这里
  const T &operator()(size_t idx, bool check = false) const noexcept {
    if (check && !isinside(idx))
      return None;
    if (bv_[idx])
      return vec_data_[rank1(idx)];
    else
      return None;
  }

  // 非 const 3D 版本 → 返回 T&（不变）
  T &operator()(size_t z, size_t y, size_t x, bool check = false) noexcept {
    if (check && !isinside(z, y, x))
      return None;
    size_t idx = ravel_idx(z, y, x);
    if (bv_[idx])
      return vec_data_[rank1(idx)];
    else
      return None;
  }

  // const 3D 版本 → 返回 const T& ← 改这里
  const T &operator()(size_t z, size_t y, size_t x,
                      bool check = false) const noexcept {
    if (check && !isinside(z, y, x))
      return None;
    size_t idx = ravel_idx(z, y, x);
    if (bv_[idx])
      return vec_data_[rank1(idx)];
    else
      return None;
  }
};