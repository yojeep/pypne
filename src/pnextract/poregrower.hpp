#include "globals.hpp"
#include "ElementGNE.hpp"
#include "types.hpp"
#include <cstddef>
#include <iostream>
#include <vector>
#include "bit_vector.hpp"

static inline const int PREFETCH_OFFSET = 32;

static inline void prefetch_voxel_neighbors_read(const int *big_array, int z,
                                                 int y, int x, int NZ, int NY,
                                                 int NX) noexcept {
  size_t base = (z * NY + y) * NX + x;
  __builtin_prefetch(&big_array[base], 0, 1);           //  M
  __builtin_prefetch(&big_array[base - 1], 0, 1);       // -X
  __builtin_prefetch(&big_array[base + 1], 0, 1);       // +X
  __builtin_prefetch(&big_array[base - NX], 0, 1);      // -Y
  __builtin_prefetch(&big_array[base + NX], 0, 1);      // +Y
  __builtin_prefetch(&big_array[base - NY * NX], 0, 1); // -Z
  __builtin_prefetch(&big_array[base + NY * NX], 0, 1); // +Z
}

static inline void prefetch_voxel_neighbors_write(const int *big_array, int z,
                                                  int y, int x, int NZ, int NY,
                                                  int NX) noexcept {
  size_t base = (z * NY + y) * NX + x;
  __builtin_prefetch(&big_array[base], 1, 1);           //  M
  __builtin_prefetch(&big_array[base - 1], 1, 1);       // -X
  __builtin_prefetch(&big_array[base + 1], 1, 1);       // +X
  __builtin_prefetch(&big_array[base - NX], 1, 1);      // -Y
  __builtin_prefetch(&big_array[base + NX], 1, 1);      // +Y
  __builtin_prefetch(&big_array[base - NY * NX], 1, 1); // -Z
  __builtin_prefetch(&big_array[base + NY * NX], 1, 1); // +Z
}

class poreNE;

enum class GrowStrategy { MedStrict, Median, MedEqs, MedEqsLoose };
class PoreGrower {
public:
  PoreGrower(TensorXXXDi32 &VElems, const BitMap3D<voxel> &vxlMap,
             const std::vector<medialBall> &ballSpace, int &firstPores,
             int &unassigned)
      : nzVE(VElems.dimension(0)), nyVE(VElems.dimension(1)),
        nxVE(VElems.dimension(2)), VElems(VElems), vxlMap(vxlMap),
        ballSpace(ballSpace), bgn(firstPores), unassigned(unassigned) {
    voxls = TensorXXXDi32(VElems.dimensions());
  }

  size_t growPores_X2() noexcept;
  size_t grow_pores() noexcept;
  void growPores() noexcept;

  template <GrowStrategy S> size_t growPoresGeneric() noexcept;
  void growPoresMedStrict() noexcept;
  void growPoresMedian() noexcept;
  void growPoresMedEqs() noexcept;
  void growPoresMedEqsLoose() noexcept;
  void retreatPoresMedian() noexcept;
  void medianElem() noexcept;
  void grow() noexcept;

private:
  const int nzVE, nyVE, nxVE;
  TensorXXXDi32 &VElems;
  TensorXXXDi32 voxls;

  const BitMap3D<voxel> &vxlMap;
  const std::vector<medialBall> &ballSpace;

  int bgn, unassigned;
};

static std::pair<int, short> get_max_count_nei(const std::array<int, 6> &neis,
                                               const int n) noexcept;

template <typename T>
inline void assign_voxel(Eigen::Tensor<T, 3, Eigen::RowMajor> &from,
                         Eigen::Tensor<T, 3, Eigen::RowMajor> &to) noexcept {
  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = from.size();

  T *from_ptr = from.data();
  T *to_ptr = to.data();

  pool.detach_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) noexcept {
        size_t count = end_idx - start_idx;
        std::memcpy(&to_ptr[start_idx], &from_ptr[start_idx],
                    count * sizeof(T));
      });

  pool.wait();
}

inline size_t PoreGrower::grow_pores() noexcept {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  size_t total_iterations = vxlMap.vec_data_.size();
  std::atomic<size_t> nChanges(0);

  pool.detach_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) noexcept {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          if (idx + PREFETCH_OFFSET < end_idx) {
            const voxel &vpre = vxlMap.vec_data_[idx + PREFETCH_OFFSET];
            prefetch_voxel_neighbors_read(voxls.data(), vpre.z + 1, vpre.y + 1,
                                          vpre.x + 1, nzVE, nyVE, nxVE);
            prefetch_voxel_neighbors_write(VElems.data(), vpre.z + 1,
                                           vpre.y + 1, vpre.x + 1, nzVE, nyVE,
                                           nxVE);
          }

          const voxel &vi = vxlMap.vec_data_[idx];
          int z = vi.z + 1, y = vi.y + 1, x = vi.x + 1;
          const int &pID = voxls(z, y, x);
          if (pID != unassigned)
            continue;

          auto try_assign = [&](auto &&v_func) noexcept -> bool {
            auto val = v_func();
            if (bgn <= val) {
              VElems(z, y, x) = val;
              ++local_nChanges;
              return true;
            }
            return false;
          };

          (void)(try_assign([&]() noexcept { return voxls(z, y, x - 1); }) ||
                 try_assign([&]() noexcept { return voxls(z, y, x + 1); }) ||
                 try_assign([&]() noexcept { return voxls(z, y - 1, x); }) ||
                 try_assign([&]() noexcept { return voxls(z, y + 1, x); }) ||
                 try_assign([&]() noexcept { return voxls(z - 1, y, x); }) ||
                 try_assign([&]() noexcept { return voxls(z + 1, y, x); }));
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });

  pool.wait();
  return nChanges.load();
}

inline size_t PoreGrower::growPores_X2() noexcept {
  size_t nChanges = 0;

  // First round
  nChanges = grow_pores();
  std::cout << "  ngrowX3:" << nChanges << ",";

  // Second round
  nChanges = grow_pores();
  std::cout << nChanges << ",";

  // Third round
  nChanges = grow_pores();
  std::cout << "  ngrowX2:" << nChanges << "  ";

  return nChanges;
}

inline void PoreGrower::growPores() noexcept {
  size_t nChanges = grow_pores();
  std::cout << "  ngrowPors:" << nChanges << "  ";
}

// =========================
// 1. RadiusMode（Compile-time）
// =========================
inline std::pair<int, short> get_max_count_nei(const std::array<int, 6> &neis,
                                               const int n) noexcept {
  int most_val = -257;
  short most_count = 0;

  for (int i = 0; i < n; ++i) {
    const int val = neis[i];

    short count = 0;
    for (int j = 0; j < n; ++j)
      count += (neis[j] == val);

    if (count > most_count || (count == most_count && val < most_val)) {
      most_count = count;
      most_val = val;
    }
  }

  return {most_val, most_count};
}

struct NeighborOffset {
  int dz, dy, dx;
};

static constexpr std::array<NeighborOffset, 6> NEIGHBOR_OFFSETS = {{
    {0, 0, -1}, // x-1
    {0, 0, 1},  // x+1
    {0, -1, 0}, // y-1
    {0, 1, 0},  // y+1
    {-1, 0, 0}, // z-1
    {1, 0, 0}   // z+1
}};

template <GrowStrategy S> size_t PoreGrower::growPoresGeneric() noexcept {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = vxlMap.vec_data_.size();
  std::atomic<size_t> nChanges(0);

  constexpr short thresh_nDiff = (S == GrowStrategy::MedStrict) ? 3 : 2;
  constexpr short thresh_most = (S == GrowStrategy::MedStrict) ? 3 : 2;

  pool.detach_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) noexcept {
        size_t local_nChanges = 0;

        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          if (idx + PREFETCH_OFFSET < end_idx) {
            const voxel &vpre = vxlMap.vec_data_[idx + PREFETCH_OFFSET];
            prefetch_voxel_neighbors_read(voxls.data(), vpre.z + 1, vpre.y + 1,
                                          vpre.x + 1, nzVE, nyVE, nxVE);
            prefetch_voxel_neighbors_write(VElems.data(), vpre.z + 1,
                                           vpre.y + 1, vpre.x + 1, nzVE, nyVE,
                                           nxVE);
          }
          const voxel &vi = vxlMap.vec_data_[idx];
          const int z = vi.z + 1, y = vi.y + 1, x = vi.x + 1;
          const float R = vi.R;
          const int &pID = voxls(z, y, x);
          if (pID != unassigned)
            continue;

          short nDifferentID = 0;
          std::array<int, 6> neis{};
          short nei_count = 0;

          for (const auto &off : NEIGHBOR_OFFSETS) {
            int zn = z + off.dz;
            int yn = y + off.dy;
            int xn = x + off.dx;

            const int &neiID = voxls(zn, yn, xn);
            if (neiID < bgn)
              continue;
            const float neiR = vxlMap(zn - 1, yn - 1, xn - 1).R;

            bool accept = false;
            if constexpr (S == GrowStrategy::MedStrict) {
              accept = (neiR >= R);
            } else if constexpr (S == GrowStrategy::Median) {
              accept = (neiR > R);
            } else if constexpr (S == GrowStrategy::MedEqs) {
              accept = (neiR >= R);
            } else if constexpr (S == GrowStrategy::MedEqsLoose) {
              accept = true;
            }

            if (!accept)
              continue;

            ++nDifferentID;

            if constexpr (S == GrowStrategy::MedStrict) {
              if (neiR > R)
                neis[nei_count++] = neiID;
            } else {
              neis[nei_count++] = neiID;
            }
          } // end neighbor loop

          if (nDifferentID >= thresh_nDiff && nei_count > 0) {
            auto [mostID, mostCNT] = get_max_count_nei(neis, nei_count);
            if (mostCNT >= thresh_most) {
              VElems(z, y, x) = mostID;
              ++local_nChanges;
            }
          }
        }

        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });
  pool.wait();
  return nChanges.load();
}

inline void PoreGrower::growPoresMedStrict() noexcept {
  size_t nChanges = PoreGrower::growPoresGeneric<GrowStrategy::MedStrict>();
  std::cout << "  ngMedStrict: " << nChanges << " ";
}

inline void PoreGrower::growPoresMedian() noexcept {
  size_t nChanges = PoreGrower::growPoresGeneric<GrowStrategy::Median>();
  std::cout << "  ngMedian: " << nChanges << " ";
}

inline void PoreGrower::growPoresMedEqs() noexcept {
  size_t nChanges = PoreGrower::growPoresGeneric<GrowStrategy::MedEqs>();
  std::cout << "  ngMedEqs: " << nChanges << " ";
}

inline void PoreGrower::growPoresMedEqsLoose() noexcept {
  size_t nChanges = PoreGrower::growPoresGeneric<GrowStrategy::MedEqsLoose>();
  std::cout << "  ngMedLoose: " << nChanges << " ";
}

inline void PoreGrower::retreatPoresMedian() noexcept {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = vxlMap.vec_data_.size();
  std::atomic<size_t> nChanges(0);

  pool.detach_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) noexcept {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          if (idx + PREFETCH_OFFSET < end_idx) {
            const voxel &vpre = vxlMap.vec_data_[idx + PREFETCH_OFFSET];
            prefetch_voxel_neighbors_read(voxls.data(), vpre.z + 1, vpre.y + 1,
                                          vpre.x + 1, nzVE, nyVE, nxVE);
            prefetch_voxel_neighbors_write(VElems.data(), vpre.z + 1,
                                           vpre.y + 1, vpre.x + 1, nzVE, nyVE,
                                           nxVE);
          }
          const voxel &vi = vxlMap.vec_data_[idx];
          const int z = vi.z + 1, y = vi.y + 1, x = vi.x + 1;
          const int &pID = voxls(z, y, x);

          if (pID < bgn)
            continue;
          short nSameID = 0;
          short nDifferentID = 0;

          auto visit_neighbor = [&](int neID) noexcept {
            if (neID == pID) {
              ++nSameID;
            } else if (neID >= bgn) {
              ++nDifferentID;
            }
          };

          visit_neighbor(voxls(z, y, x - 1));
          visit_neighbor(voxls(z, y, x + 1));
          visit_neighbor(voxls(z, y - 1, x));
          visit_neighbor(voxls(z, y + 1, x));
          visit_neighbor(voxls(z - 1, y, x));
          visit_neighbor(voxls(z + 1, y, x));

          if (nDifferentID > 0 && nSameID > 0) {
            VElems(z, y, x) = unassigned;
            ++local_nChanges;
          }
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });
  pool.wait();

  (std::cout << "  nRetreat:" << nChanges.load()).flush();
}

inline void PoreGrower::medianElem() noexcept {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = vxlMap.vec_data_.size();
  std::atomic<size_t> nChanges{0};

  pool.detach_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) noexcept {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          if (idx + PREFETCH_OFFSET < end_idx) {
            const voxel &vpre = vxlMap.vec_data_[idx + PREFETCH_OFFSET];
            prefetch_voxel_neighbors_read(voxls.data(), vpre.z + 1, vpre.y + 1,
                                          vpre.x + 1, nzVE, nyVE, nxVE);
            prefetch_voxel_neighbors_write(VElems.data(), vpre.z + 1,
                                           vpre.y + 1, vpre.x + 1, nzVE, nyVE,
                                           nxVE);
          }
          const voxel &vi = vxlMap.vec_data_[idx];
          const int z = vi.z + 1, y = vi.y + 1, x = vi.x + 1;
          const int &pID = voxls(z, y, x);

          if (pID < bgn)
            continue;

          short nSameID = 0;
          short nDifferentID = 0;
          std::array<int, 6> neis{};
          short nei_count = 0;

          auto visit_neighbor = [&](int neID) noexcept {
            if (neID == pID) {
              ++nSameID;
            } else if (neID >= bgn) {
              ++nDifferentID;
              neis[nei_count++] = neID;
            }
          };

          visit_neighbor(voxls(z, y, x - 1));
          visit_neighbor(voxls(z, y, x + 1));
          visit_neighbor(voxls(z, y - 1, x));
          visit_neighbor(voxls(z, y + 1, x));
          visit_neighbor(voxls(z - 1, y, x));
          visit_neighbor(voxls(z + 1, y, x));

          if (nDifferentID <= nSameID || nei_count == 0)
            continue;

          auto [best_id, best_cnt] = get_max_count_nei(neis, nei_count);

          if (best_cnt > nSameID) {
            VElems(z, y, x) = best_id;
            ++local_nChanges;
          }
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });
  pool.wait();
  std::cout << "  nMedian:" << nChanges.load() << " ";
}

inline void PoreGrower::grow() noexcept {
  growPoresMedStrict();
  growPoresMedStrict();
  growPoresMedStrict();
  growPoresMedian();
  growPoresMedStrict();
  growPoresMedian();
  growPoresMedStrict();
  growPoresMedian();
  growPoresMedStrict();
  std::cout << std::endl;
  growPoresMedian();
  growPoresMedStrict();
  growPoresMedian();
  growPoresMedStrict();
  growPoresMedian();
  growPoresMedStrict();
  std::cout << std::endl;
  growPoresMedEqs();
  growPoresMedian();
  growPoresMedStrict();
  growPoresMedian();
  growPoresMedian();
  growPoresMedStrict();
  growPoresMedEqs();
  growPoresMedian();
  growPoresMedian();
  growPoresMedStrict();
  std::cout << std::endl;
  growPoresMedEqs();
  growPoresMedEqs();
  growPoresMedEqs();
  growPoresMedEqs();
  growPoresMedEqs();
  growPoresMedEqsLoose();
  growPoresMedEqs();
  std::cout << std::endl;
  growPoresMedEqs();
  growPoresMedEqsLoose();
  growPoresMedEqs();
  growPoresMedEqsLoose();
  growPoresMedEqs();
  growPoresMedEqsLoose();
  growPoresMedEqs();
  growPoresMedEqsLoose();
  growPoresMedEqs();
  growPoresMedEqsLoose();
  growPoresMedEqs();
  growPoresMedEqsLoose();

  std::cout << std::endl;

  growPores();
  growPoresMedian();
  growPoresMedEqs();
  growPores();
  growPoresMedian();
  growPores();
  growPores();
  growPores();
  std::cout << std::endl;
  growPores();
  growPores();
  growPores();
  growPoresMedian();
  growPores();
  growPores();
  growPores();
  std::cout << std::endl;
  growPores();
  growPores();
  growPores();
  growPores();
  growPores();
  growPores();
  std::cout << std::endl;
  growPores();
  growPores_X2();
  growPores_X2();
  growPores_X2();
  growPores_X2();
  std::cout << std::endl;

  medianElem();
  medianElem();
  medianElem();

  medianElem();
  medianElem();

  std::cout << std::endl;

  growPores();
  while (growPores_X2())
    ;
  growPores();

  std::cout << std::endl;

  retreatPoresMedian();

  for (const medialBall &bi : ballSpace) {
    const medialBall *mastrSphere = bi.mastrSphere();
    int izi = bi.iz() + 1;
    int iyi = bi.iy() + 1;
    int ixi = bi.ix() + 1;

    int izm = mastrSphere->iz() + 1;
    int iym = mastrSphere->iy() + 1;
    int ixm = mastrSphere->ix() + 1;
    VElems(izi, iyi, ixi) = VElems(izm, iym, ixm);
  }
  growPoresMedian();
  growPoresMedEqs();
  growPoresMedEqs();
  growPoresMedEqs();
  growPoresMedian();
  growPoresMedEqs();
  growPoresMedEqs();
  growPoresMedEqsLoose();
  growPoresMedEqs();
  growPoresMedEqsLoose();
  growPoresMedEqs();
  growPoresMedEqsLoose();
  growPoresMedEqs();
  growPores();
  growPores_X2();
  growPoresMedEqs();
  growPores();
  growPores_X2();

  std::cout << std::endl;

  medianElem();
  for (const medialBall &bi : ballSpace) {
    int izi = bi.iz() + 1;
    int iyi = bi.iy() + 1;
    int ixi = bi.ix() + 1;
    const medialBall *mastrSphere = bi.mastrSphere();
    int izm = mastrSphere->iz() + 1;
    int iym = mastrSphere->iy() + 1;
    int ixm = mastrSphere->ix() + 1;

    VElems(izi, iyi, ixi) = VElems(izm, iym, ixm);
  }
  growPoresMedEqsLoose();
  int label = bgn;
  for (const auto &bi : ballSpace) {
    if (bi.boss == &bi) {
      int izi = bi.iz() + 1;
      int iyi = bi.iy() + 1;
      int ixi = bi.ix() + 1;
      VElems(izi, iyi, ixi) = label;
      ++label;
    }
  }

  std::cout << std::endl;
}