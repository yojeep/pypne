#include "blockNet.h"
#include <atomic>
#include <cassert>
#include <cstddef>
#include <vector>

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

size_t grow_pores(voxelField<int> &VElems, voxelField<int> &voxls, int &bgn,
                  int &lst, int &unassigned) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny(), nx = VElems.nx();
  const size_t total_iterations = static_cast<size_t>(nz * ny * nx);
  std::atomic<size_t> nChanges(0);
  IndexUnraveler unraveler({static_cast<size_t>(nz), static_cast<size_t>(ny),
                            static_cast<size_t>(nx)});

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        size_t local_nChanges = 0;
        std::array<size_t, 3> coords{};
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          unraveler.unravel(idx, coords);
          int z = coords[0], y = coords[1], x = coords[2];
          if (z == 0 || z == nz - 1 || y == 0 || y == ny - 1 || x == 0 ||
              x == nx - 1)
            continue;

          const int *pijk = &voxls(idx);
          if (*pijk != unassigned)
            continue;

          auto try_assign = [&](auto &&v_func) -> bool {
            auto val = v_func(pijk);
            if (bgn <= val) {
              VElems(idx) = val;
              ++local_nChanges;
              return true;
            }
            return false;
          };

          (void)(try_assign([&](const auto &p) { return voxls.v_i(1, p); }) ||
                 try_assign([&](const auto &p) { return voxls.v_i(-1, p); }) ||
                 try_assign([&](const auto &p) { return voxls.v_j(1, p); }) ||
                 try_assign([&](const auto &p) { return voxls.v_j(-1, p); }) ||
                 try_assign([&](const auto &p) { return voxls.v_k(1, p); }) ||
                 try_assign([&](const auto &p) { return voxls.v_k(-1, p); }));
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });

  pool.wait();
  return nChanges.load();
}

size_t growPores_X2(voxelField<int> &VElems, voxelField<int> &voxls, int &bgn,
                    int &lst, int &unassigned) {
  size_t nChanges = 0;

  // First round
  nChanges = grow_pores(VElems, voxls, bgn, lst, unassigned);
  std::cout << "  ngrowX3:" << nChanges << ",";

  // Second round
  nChanges = grow_pores(VElems, voxls, bgn, lst, unassigned);
  std::cout << nChanges << ",";

  // Third round
  nChanges = grow_pores(VElems, voxls, bgn, lst, unassigned);
  std::cout << "  ngrowX2:" << nChanges << "  ";

  return nChanges;
}

void growPores(voxelField<int> &VElems, voxelField<int> &voxls, int &bgn,
               int &lst, int &unassigned) {
  size_t nChanges = grow_pores(VElems, voxls, bgn, lst, unassigned);
  std::cout << "  ngrowPors:" << nChanges << "  ";
}

// =========================
// 1. RadiusMode（Compile-time）
// =========================
inline std::pair<int, short> get_max_count_nei(const std::array<int, 6> &neis,
                                               int n) noexcept {
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

enum class GrowStrategy { MedStrict, Median, MedEqs, MedEqsLoose };

struct NeighborOffset {
  int dk, dj, di;
};

static constexpr std::array<NeighborOffset, 6> NEIGHBOR_OFFSETS = {{
    {0, 0, -1}, // i-1
    {0, 0, 1},  // i+1
    {0, -1, 0}, // j-1
    {0, 1, 0},  // j+1
    {-1, 0, 0}, // k-1
    {1, 0, 0}   // k+1
}};

template <GrowStrategy S>
size_t growPoresGeneric(const inputDataNE &cg, voxelField<int> &VElems,
                        voxelField<int> &voxls, int bgn, int lst,
                        const std::vector<poreNE *> &poreIs, int unassigned) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny();
  const size_t total_iterations = static_cast<size_t>(nz) * ny;
  std::atomic<size_t> nChanges(0);
  IndexUnraveler unraveler({static_cast<size_t>(nz), static_cast<size_t>(ny)});

  constexpr short thresh_nDiff = (S == GrowStrategy::MedStrict) ? 3 : 2;
  constexpr short thresh_most = (S == GrowStrategy::MedStrict) ? 3 : 2;

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        size_t local_nChanges = 0;
        std::array<size_t, 2> coords{};

        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          unraveler.unravel(idx, coords);
          const int k = coords[0], j = coords[1];
          if (k == 0 || k == nz - 1 || j == 0 || j == ny - 1)
            continue;

          const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
          for (int ix = 0; ix < s.cnt; ++ix) {
            for (int i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i) {
              const int *pijk = &voxls(i, j, k);
              if (*pijk != unassigned)
                continue;

              const float R = cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 1)->R;
              short nDifferentID = 0;
              std::array<int, 6> neis{};
              short nei_count = 0;

              for (const auto &off : NEIGHBOR_OFFSETS) {
                const int nk = k + off.dk;
                const int nj = j + off.dj;
                const int ni = i + off.di;

                const int neiID = voxls(ni, nj, nk);
                if (neiID < bgn)
                  continue;

                const float neiR =
                    cg.segs_[(nk - 1) * cg.ny + (nj - 1)].vxl(ni - 1)->R;

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
                  VElems(i, j, k) = mostID;
                  ++local_nChanges;
                }
              }
            }
          }
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });
  pool.wait();
  return nChanges.load();
}

void growPoresMedStrict(const inputDataNE &cg, voxelField<int> &VElems,
                        voxelField<int> &voxls, int &bgn, int &lst,
                        const std::vector<poreNE *> &poreIs, int &unassigned) {
  size_t nChanges = growPoresGeneric<GrowStrategy::MedStrict>(
      cg, VElems, voxls, bgn, lst, poreIs, unassigned);
  std::cout << "  ngMedStrict: " << nChanges << " ";
}

void growPoresMedian(const inputDataNE &cg, voxelField<int> &VElems,
                     voxelField<int> &voxls, int &bgn, int &lst,
                     const std::vector<poreNE *> &poreIs, int &unassigned) {
  size_t nChanges = growPoresGeneric<GrowStrategy::Median>(
      cg, VElems, voxls, bgn, lst, poreIs, unassigned);
  std::cout << "  ngMedian: " << nChanges << " ";
}

void growPoresMedEqs(const inputDataNE &cg, voxelField<int> &VElems,
                     voxelField<int> &voxls, int &bgn, int &lst,
                     const std::vector<poreNE *> &poreIs, int &unassigned) {
  size_t nChanges = growPoresGeneric<GrowStrategy::MedEqs>(
      cg, VElems, voxls, bgn, lst, poreIs, unassigned);
  std::cout << "  ngMedEqs: " << nChanges << " ";
}

void growPoresMedEqsLoose(const inputDataNE &cg, voxelField<int> &VElems,
                          voxelField<int> &voxls, int &bgn, int &lst,
                          const std::vector<poreNE *> &poreIs,
                          int &unassigned) {
  size_t nChanges = growPoresGeneric<GrowStrategy::MedEqsLoose>(
      cg, VElems, voxls, bgn, lst, poreIs, unassigned);
  std::cout << "  ngMedLoose: " << nChanges << " ";
}

void retreatPoresMedian(const inputDataNE &cg, voxelField<int> &VElems,
                        voxelField<int> &voxls, int &bgn, int &lst,
                        const std::vector<poreNE *> &poreIs, int &unassigned) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny();
  const size_t total_iterations = static_cast<size_t>(nz * ny);
  std::atomic<size_t> nChanges(0);
  IndexUnraveler unraveler({static_cast<size_t>(nz), static_cast<size_t>(ny)});

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        size_t local_nChanges = 0;
        std::array<size_t, 2> coords{};
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          unraveler.unravel(idx, coords);
          int k = coords[0], j = coords[1];
          if (k == 0 || k == nz - 1 || j == 0 || j == ny - 1)
            continue;
          const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
          for (int ix = 0; ix < s.cnt; ++ix) {
            for (int i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i) {
              const int *pijk = &voxls(i, j, k);

              int pID = *pijk;
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

              visit_neighbor(voxls.v_i(-1, pijk));
              visit_neighbor(voxls.v_i(1, pijk));
              visit_neighbor(voxls.v_j(-1, pijk));
              visit_neighbor(voxls.v_j(1, pijk));
              visit_neighbor(voxls.v_k(-1, pijk));
              visit_neighbor(voxls.v_k(1, pijk));

              if (nDifferentID > 0 && nSameID > 0) {
                VElems(i, j, k) = unassigned;
                ++local_nChanges;
              }
            }
          }
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });
  pool.wait();

  (std::cout << "  nRetreat:" << nChanges.load()).flush();
}

void medianElem(const inputDataNE &cg, voxelField<int> &VElems,
                voxelField<int> &voxls, int &bgn, int &lst,
                const std::vector<poreNE *> &poreIs) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny();
  const size_t total_iterations = static_cast<size_t>(nz * ny);
  IndexUnraveler unraveler({static_cast<size_t>(nz), static_cast<size_t>(ny)});
  std::atomic<size_t> nChanges{0};

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        std::array<size_t, 2> coords{};
        size_t local_nChanges = 0;

        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          unraveler.unravel(idx, coords);
          const int k = coords[0], j = coords[1];
          if (k == 0 || k == nz - 1 || j == 0 || j == ny - 1)
            continue;
          const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
          for (int ix = 0; ix < s.cnt; ++ix) {
            for (int i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i) {
              const int *pijk = &voxls(i, j, k);
              const long pID = *pijk;

              if (pID < bgn)
                continue;

              short nSameID = 0;
              short nDifferentID = 0;
              std::array<int, 6> nei_buf{};
              short nei_count = 0;

              auto visit_neighbor = [&](long neID) {
                if (neID == pID) {
                  ++nSameID;
                } else if (neID >= bgn) {
                  ++nDifferentID;
                  nei_buf[nei_count++] = static_cast<int>(neID);
                }
              };

              visit_neighbor(voxls.v_i(-1, pijk));
              visit_neighbor(voxls.v_i(1, pijk));
              visit_neighbor(voxls.v_j(-1, pijk));
              visit_neighbor(voxls.v_j(1, pijk));
              visit_neighbor(voxls.v_k(-1, pijk));
              visit_neighbor(voxls.v_k(1, pijk));

              if (nDifferentID <= nSameID || nei_count == 0)
                continue;

              auto [best_id, best_cnt] = get_max_count_nei(nei_buf, nei_count);

              if (best_cnt > nSameID) {
                ++local_nChanges;
                VElems(i, j, k) = best_id;
              }
            }
          }
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });
  pool.wait();
  std::cout << "  nMedian:" << nChanges.load() << " ";
}
