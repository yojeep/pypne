#include "Eigen/CXX11/Tensor"
#include "ElementGNE.h"
#include "inputData.h"
#include "voxelImage.h"
#include <cstddef>
#include <iostream>
#include <vector>

// bool parseBool(std::string_view sv) {
//     if (sv.empty()) return false;
//     char c = sv[0];
//     return (c == 'T' || c == 't' || c == 'Y' || c == 'y' || c == '1');
// }

class poreNE;

enum class GrowStrategy { MedStrict, Median, MedEqs, MedEqsLoose };
class PoreGrower {
public:
  PoreGrower(const inputDataNE &cg, const std::vector<poreNE *> &poreIs,
             voxelField<int> &VElems, voxelField<int> &voxls,
             const std::vector<voxel> &vxlSpace,
             const Eigen::Tensor<voxel *, 3, Eigen::RowMajor> &vxlMap,
             const std::vector<medialBall> &ballSpace, int &firstPores,
             int &lastPores, int &unassigned)
      : cg(cg), poreIs(poreIs), VElems(VElems), voxls(voxls),
        vxlSpace(vxlSpace), vxlMap(vxlMap), ballSpace(ballSpace),
        bgn(firstPores), lst(lastPores), unassigned(unassigned) {}

  size_t growPores_X2();
  size_t grow_pores();
  void growPores();

  template <GrowStrategy S> size_t growPoresGeneric();
  void growPoresMedStrict();
  void growPoresMedian();
  void growPoresMedEqs();
  void growPoresMedEqsLoose();
  void retreatPoresMedian();
  void medianElem();
  void grow();
  //   inline std::pair<int, short> get_max_count_nei(const std::array<int, 6>
  //   &neis,
  //                                                  const int n) noexcept;
  //   int unassigned() const { return unassigned_; }

private:
  const inputDataNE &cg;
  const std::vector<poreNE *> &poreIs;
  voxelField<int> &VElems;
  voxelField<int> &voxls;
  const std::vector<voxel> &vxlSpace;
  const Eigen::Tensor<voxel *, 3, Eigen::RowMajor> &vxlMap;
  const std::vector<medialBall> &ballSpace;

  int &bgn, &lst, &unassigned;
};

static std::pair<int, short> get_max_count_nei(const std::array<int, 6> &neis,
                                               const int n) noexcept;

template <typename T>
inline void assign_voxel(voxelField<T> &from, voxelField<T> &to) {
  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = from.size();

  T *from_ptr = from.data_.data();
  T *to_ptr = to.data_.data();

  pool.detach_blocks(0, total_iterations,
                     [&](const size_t start_idx, const size_t end_idx) {
                       size_t count = end_idx - start_idx;
                       std::memcpy(&to_ptr[start_idx], &from_ptr[start_idx],
                                   count * sizeof(T));
                     });

  pool.wait();
}

inline size_t PoreGrower::grow_pores() {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  size_t total_iterations = vxlSpace.size();
  std::atomic<size_t> nChanges(0);

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          const voxel &vi = vxlSpace[idx];
          int z = vi.k + 1, y = vi.j + 1, x = vi.i + 1;
          const int &pID = voxls(x, y, z);
          if (pID != unassigned)
            continue;

          auto try_assign = [&](auto &&v_func) -> bool {
            auto val = v_func(pID);
            if (bgn <= val) {
              VElems(x, y, z) = val;
              ++local_nChanges;
              return true;
            }
            return false;
          };

          (void)(try_assign(
                     [&](const auto &p) { return voxls(x - 1, y, z); }) ||
                 try_assign(
                     [&](const auto &p) { return voxls(x + 1, y, z); }) ||
                 try_assign(
                     [&](const auto &p) { return voxls(x, y - 1, z); }) ||
                 try_assign(
                     [&](const auto &p) { return voxls(x, y + 1, z); }) ||
                 try_assign(
                     [&](const auto &p) { return voxls(x, y, z - 1); }) ||
                 try_assign([&](const auto &p) { return voxls(x, y, z + 1); }));
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });

  pool.wait();
  return nChanges.load();
}

inline size_t PoreGrower::growPores_X2() {
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

inline void PoreGrower::growPores() {
  size_t nChanges = grow_pores();
  std::cout << "  ngrowPors:" << nChanges << "  ";
}

// =========================
// 1. RadiusMode（Compile-time）
// =========================
inline std::pair<int, short> get_max_count_nei(const std::array<int, 6> &neis,
                                               const int n) noexcept {
  int most_val = -1;
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
    {0, 0, -1}, // i-1
    {0, 0, 1},  // i+1
    {0, -1, 0}, // j-1
    {0, 1, 0},  // j+1
    {-1, 0, 0}, // k-1
    {1, 0, 0}   // k+1
}};

template <GrowStrategy S> size_t PoreGrower::growPoresGeneric() {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = vxlSpace.size();
  std::atomic<size_t> nChanges(0);

  constexpr short thresh_nDiff = (S == GrowStrategy::MedStrict) ? 3 : 2;
  constexpr short thresh_most = (S == GrowStrategy::MedStrict) ? 3 : 2;

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        size_t local_nChanges = 0;

        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          const voxel &vi = vxlSpace[idx];
          const int z = vi.k + 1, y = vi.j + 1, x = vi.i + 1;
          const float R = vi.R;
          const int &pID = voxls(x, y, z);
          if (pID != unassigned)
            continue;

          short nDifferentID = 0;
          std::array<int, 6> neis{};
          short nei_count = 0;

          for (const auto &off : NEIGHBOR_OFFSETS) {
            int zn = z + off.dz;
            int yn = y + off.dy;
            int xn = x + off.dx;

            const int &neiID = voxls(xn, yn, zn);
            if (neiID < bgn)
              continue;
            const float neiR = vxlMap(zn - 1, yn - 1, xn - 1)->R;

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
              VElems(x, y, z) = mostID;
              ++local_nChanges;
            }
          }
        }

        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });
  pool.wait();
  return nChanges.load();
}

inline void PoreGrower::growPoresMedStrict() {
  size_t nChanges = PoreGrower::growPoresGeneric<GrowStrategy::MedStrict>();
  std::cout << "  ngMedStrict: " << nChanges << " ";
}

inline void PoreGrower::growPoresMedian() {
  size_t nChanges = PoreGrower::growPoresGeneric<GrowStrategy::Median>();
  std::cout << "  ngMedian: " << nChanges << " ";
}

inline void PoreGrower::growPoresMedEqs() {
  size_t nChanges = PoreGrower::growPoresGeneric<GrowStrategy::MedEqs>();
  std::cout << "  ngMedEqs: " << nChanges << " ";
}

inline void PoreGrower::growPoresMedEqsLoose() {
  size_t nChanges = PoreGrower::growPoresGeneric<GrowStrategy::MedEqsLoose>();
  std::cout << "  ngMedLoose: " << nChanges << " ";
}

inline void PoreGrower::retreatPoresMedian() {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = vxlSpace.size();
  std::atomic<size_t> nChanges(0);

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          const voxel &vi = vxlSpace[idx];
          const int z = vi.k + 1, y = vi.j + 1, x = vi.i + 1;
          const int &pID = voxls(x, y, z);

          if (pID < bgn)
            continue;
          short nSameID = 0;
          short nDifferentID = 0;

          auto visit_neighbor = [&](int neID) {
            if (neID == pID) {
              ++nSameID;
            } else if (neID >= bgn) {
              ++nDifferentID;
            }
          };

          visit_neighbor(voxls(x - 1, y, z));
          visit_neighbor(voxls(x + 1, y, z));
          visit_neighbor(voxls(x, y - 1, z));
          visit_neighbor(voxls(x, y + 1, z));
          visit_neighbor(voxls(x, y, z - 1));
          visit_neighbor(voxls(x, y, z + 1));

          if (nDifferentID > 0 && nSameID > 0) {
            VElems(x, y, z) = unassigned;
            ++local_nChanges;
          }
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });
  pool.wait();

  (std::cout << "  nRetreat:" << nChanges.load()).flush();
}

inline void PoreGrower::medianElem() {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = vxlSpace.size();
  std::atomic<size_t> nChanges{0};

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          const voxel &vi = vxlSpace[idx];
          const int z = vi.k + 1, y = vi.j + 1, x = vi.i + 1;
          const int &pID = voxls(x, y, z);

          if (pID < bgn)
            continue;

          short nSameID = 0;
          short nDifferentID = 0;
          std::array<int, 6> neis{};
          short nei_count = 0;

          auto visit_neighbor = [&](long neID) {
            if (neID == pID) {
              ++nSameID;
            } else if (neID >= bgn) {
              ++nDifferentID;
              neis[nei_count++] = static_cast<int>(neID);
            }
          };

          visit_neighbor(voxls(x - 1, y, z));
          visit_neighbor(voxls(x + 1, y, z));
          visit_neighbor(voxls(x, y - 1, z));
          visit_neighbor(voxls(x, y + 1, z));
          visit_neighbor(voxls(x, y, z - 1));
          visit_neighbor(voxls(x, y, z + 1));

          if (nDifferentID <= nSameID || nei_count == 0)
            continue;

          auto [best_id, best_cnt] = get_max_count_nei(neis, nei_count);

          if (best_cnt > nSameID) {
            ++local_nChanges;
            VElems(x, y, z) = best_id;
          }
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });
  pool.wait();
  std::cout << "  nMedian:" << nChanges.load() << " ";
}

inline void PoreGrower::grow() {
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

  for (const auto &bi : ballSpace) {
    medialBall *mastrSphere = bi.mastrSphere();
    VElems(bi.fi + 1, bi.fj + 1, bi.fk + 1) =
        VElems(mastrSphere->fi + 1, mastrSphere->fj + 1, mastrSphere->fk + 1);
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
  for (const auto &bi : ballSpace) {
    medialBall *mastrSphere = bi.mastrSphere();
    VElems(bi.fi + 1, bi.fj + 1, bi.fk + 1) =
        VElems(mastrSphere->fi + 1, mastrSphere->fj + 1, mastrSphere->fk + 1);
  }
  growPoresMedEqsLoose();
  int label = bgn;
  for (const auto &bi : ballSpace) {
    if (bi.boss == &bi) {
      VElems(bi.fi + 1, bi.fj + 1, bi.fk + 1) = label;
      ++label;
    }
  }

  std::cout << std::endl;
}