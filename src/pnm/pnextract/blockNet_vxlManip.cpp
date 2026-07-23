#ifndef voxelImageManip_H
#define voxelImageManip_H
#include "blockNet.h"
#include <atomic>
#include <cassert>
#include <vector>

// #include "global.h"
// class mapComparer  {  public:
// bool operator() (pair<const int,short>& i1, pair<const int,short> i2) {return
// i1.second<i2.second;}  };

using namespace std; // std::pair, vector map

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

size_t grow_pores(voxelField<int> &VElems, voxelField<int> &voxls, int bgn,
                  int lst, int porValue) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny(), nx = VElems.nx();
  const size_t total_iterations = static_cast<size_t>(nz * ny * nx);
  std::atomic<size_t> nChanges(0);
  IndexUnraveler unraveler({static_cast<size_t>(nz), static_cast<size_t>(ny),
                            static_cast<size_t>(nx)});

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        std::array<size_t, 3> coords{};
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          unraveler.unravel(idx, coords);
          int z = coords[0], y = coords[1], x = coords[2];
          if (z == 0 || z == nz - 1 || y == 0 || y == ny - 1 || x == 0 ||
              x == nx - 1)
            continue;

          const int *pijk = &voxls(idx);
          if (*pijk == porValue) {
            auto try_assign = [&](auto &&v_func) -> bool {
              auto val = v_func(pijk);
              if (bgn <= val && val <= lst) {
                VElems(idx) = val;
                ++local_nChanges;
                return true;
              }
              return false;
            };

            (void)(try_assign([&](const auto &p) { return voxls.v_i(1, p); }) ||
                   try_assign(
                       [&](const auto &p) { return voxls.v_i(-1, p); }) ||
                   try_assign([&](const auto &p) { return voxls.v_j(1, p); }) ||
                   try_assign(
                       [&](const auto &p) { return voxls.v_j(-1, p); }) ||
                   try_assign([&](const auto &p) { return voxls.v_k(1, p); }) ||
                   try_assign([&](const auto &p) { return voxls.v_k(-1, p); }));
          }
        }
        nChanges.fetch_add(local_nChanges, std::memory_order_relaxed);
      });

  pool.wait();
  return nChanges.load();
}

size_t growPores_X2(voxelField<int> &VElems, voxelField<int> &voxls, int bgn,
                    int lst, int porValue) {
  size_t nChanges = 0;

  // First round
  nChanges = grow_pores(VElems, voxls, bgn, lst, porValue);
  std::cout << "  ngrowX3:" << nChanges << ",";

  // Second round
  nChanges = grow_pores(VElems, voxls, bgn, lst, porValue);
  std::cout << nChanges << ",";

  // Third round
  nChanges = grow_pores(VElems, voxls, bgn, lst, porValue);
  std::cout << "  ngrowX2:" << nChanges << "  ";

  return nChanges;
}

void growPores(voxelField<int> &VElems, voxelField<int> &voxls, int bgn,
               int lst, int porValue) {
  size_t nChanges = grow_pores(VElems, voxls, bgn, lst, porValue);
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
void retreatPoresMedian(const inputDataNE &cg, voxelField<int> &VElems,
                        voxelField<int> &voxls, long bgn, long lst,
                        const vector<poreNE *> &poreIs, long unassigned) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny();
  // 内部区域
  const int j_start = 1, j_end = ny - 1;
  const int k_start = 1, k_end = nz - 1;
  const int j_size = j_end - j_start;
  const int k_size = k_end - k_start;
  const size_t total_kj = static_cast<size_t>(k_size * j_size); // collapse(2)
  auto futures = pool.submit_blocks(
      0, total_kj, [&](const size_t start_idx, const size_t end_idx) -> size_t {
        size_t local_nChanges = 0;

        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          // 一维 idx -> 二维 (k, j)
          int j = j_start + (idx % j_size);
          int k = k_start + (idx / j_size);
          const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
          for (short ix = 0; ix < s.cnt; ++ix) {
            for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i) {
              const int *pijk = &voxls(i, j, k);
              long pID = *pijk;
              short nSameID = 0;
              short nDifferentID = 0;

              if (pID >= bgn && pID <= lst) {
                if (voxls.v_i(-1, pijk) == pID)
                  nSameID++;
                else if (bgn <= voxls.v_i(-1, pijk) &&
                         lst >= voxls.v_i(-1, pijk))
                  nDifferentID++;
                if (voxls.v_i(1, pijk) == pID)
                  nSameID++;
                else if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk))
                  nDifferentID++;
                if (voxls.v_j(-1, pijk) == pID)
                  nSameID++;
                else if (bgn <= voxls.v_j(-1, pijk) &&
                         lst >= voxls.v_j(-1, pijk))
                  nDifferentID++;
                if (voxls.v_j(1, pijk) == pID)
                  nSameID++;
                else if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk))
                  nDifferentID++;
                if (voxls.v_k(-1, pijk) == pID)
                  nSameID++;
                else if (bgn <= voxls.v_k(-1, pijk) &&
                         lst >= voxls.v_k(-1, pijk))
                  nDifferentID++;
                if (voxls.v_k(1, pijk) == pID)
                  nSameID++;
                else if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk))
                  nDifferentID++;

                if (nDifferentID > 0 && nSameID > 0) {
                  VElems(i, j, k) = unassigned;
                  ++local_nChanges;
                }
              }
            }
          }
        }
        return local_nChanges;
      });
  size_t nChanges(0);
  for (auto &future : futures) {
    nChanges += future.get();
  }
  (cout << "  nRetreat:" << nChanges).flush();
}

void growPoresMedStrict(const inputDataNE &cg, voxelField<int> &VElems,
                        voxelField<int> &voxls, long bgn, long lst,
                        const vector<poreNE *> &poreIs, long rawValue) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny();
  const int j_start = 1, j_end = ny - 1;
  const int k_start = 1, k_end = nz - 1;
  const int j_size = j_end - j_start;
  const int k_size = k_end - k_start;
  const size_t total_iterations = static_cast<size_t>(k_size * j_size);
  auto nChanges_futures = pool.submit_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) -> size_t {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          int j = j_start + (idx % j_size);
          int k = k_start + (idx / j_size);
          const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
          for (short ix = 0; ix < s.cnt; ++ix) {
            for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i) {
              long pID = voxls(i, j, k);
              const int *pijk = &voxls(i, j, k);

              if (pID == rawValue) {
                float R = cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 1)->R;
                short nDifferentID = 0;
                if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk) &&
                    cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk) &&
                    cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
                  nDifferentID++;

                if (nDifferentID >= 3) {
                  map<int, short> neis;

                  long neI = voxls.v_i(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_i(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_j(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_j(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_k(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_k(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);

                  map<int, short>::iterator neitr =
                      max_element(neis.begin(), neis.end(), mapComparer<int>());
                  if (neitr->second >= 3) {
                    ++local_nChanges;
                    VElems(i, j, k) = neitr->first;
                  }
                }
              }
            }
          }
        }
        return local_nChanges;
      });
  size_t nChanges(0);
  for (auto &future : nChanges_futures) {
    nChanges += future.get();
  }

  cout << "  ngMedStrict:" << nChanges << " ";
}

void growPoresMedian(const inputDataNE &cg, voxelField<int> &VElems,
                     voxelField<int> &voxls, long bgn, long lst,
                     const vector<poreNE *> &poreIs, long rawValue) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny();
  const int j_start = 1, j_end = ny - 1;
  const int k_start = 1, k_end = nz - 1;
  const int j_size = j_end - j_start;
  const int k_size = k_end - k_start;
  const size_t total_iterations = static_cast<size_t>(k_size * j_size);
  auto nChanges_futures = pool.submit_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) -> size_t {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          int j = j_start + (idx % j_size);
          int k = k_start + (idx / j_size);
          const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
          for (short ix = 0; ix < s.cnt; ++ix) {
            for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i) {
              // voxel* v=s.s[ix].v(i-1);
              // if (v)

              const int *pijk = &voxls(i, j, k);
              long pID = *pijk;

              if (pID == rawValue) {
                float R = cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 1)->R;

                short nDifferentID = 0;
                if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R > R)
                  nDifferentID++;
                if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R > R)
                  nDifferentID++;
                if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R > R)
                  nDifferentID++;
                if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R > R)
                  nDifferentID++;
                if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk) &&
                    cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R > R)
                  nDifferentID++;
                if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk) &&
                    cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R > R)
                  nDifferentID++;

                if (nDifferentID >= 2) {
                  map<int, short> neis;

                  long neI = voxls.v_i(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_i(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_j(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_j(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_k(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_k(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R > R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);

                  map<int, short>::iterator neitr =
                      max_element(neis.begin(), neis.end(), mapComparer<int>());
                  if (neitr->second >= 2) {
                    ++local_nChanges;
                    VElems(i, j, k) = neitr->first;
                  }
                }
              }
            }
          }
        }
        return local_nChanges;
      });
  size_t nChanges(0);
  for (auto &future : nChanges_futures) {
    nChanges += future.get();
  }
  cout << "  ngMedian:" << nChanges << " ";
}

void growPoresMedEqs(const inputDataNE &cg, voxelField<int> &VElems,
                     voxelField<int> &voxls, long bgn, long lst,
                     const vector<poreNE *> &poreIs, long rawValue) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny();
  const int j_start = 1, j_end = ny - 1;
  const int k_start = 1, k_end = nz - 1;
  const int j_size = j_end - j_start;
  const int k_size = k_end - k_start;
  const size_t total_iterations = static_cast<size_t>(k_size * j_size);
  auto nChanges_futures = pool.submit_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) -> size_t {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          int j = j_start + (idx % j_size);
          int k = k_start + (idx / j_size);
          const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
          for (short ix = 0; ix < s.cnt; ++ix) {
            for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i) {
              // voxel* v=s.s[ix].v(i-1);
              // if (v)

              const int *pijk = &voxls(i, j, k);
              long pID = *pijk;

              if (pID == rawValue) {
                float R = cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 1)->R;

                short nDifferentID = 0;
                if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk) &&
                    cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk) &&
                    cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
                  nDifferentID++;
                if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk) &&
                    cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
                  nDifferentID++;

                if (nDifferentID >= 2) {
                  map<int, short> neis;

                  long neI = voxls.v_i(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R >= R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_i(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R >= R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_j(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R >= R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_j(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R >= R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_k(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_k(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI &&
                      cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);

                  map<int, short>::iterator neitr =
                      max_element(neis.begin(), neis.end(), mapComparer<int>());
                  if (neitr->second >= 2) {
                    ++local_nChanges;
                    VElems(i, j, k) = neitr->first;
                  }
                }
              }
            }
          }
        }
        return local_nChanges;
      });
  size_t nChanges(0);
  for (auto &future : nChanges_futures) {
    nChanges += future.get();
  }

  cout << "  ngMedEqs:" << nChanges << "  ";
}

void growPoresMedEqsLoose(const inputDataNE &cg, voxelField<int> &VElems,
                          voxelField<int> &voxls, long bgn, long lst,
                          const vector<poreNE *> &poreIs, long rawValue) {
  assign_voxel(VElems, voxls);
  auto &pool = GlobalThreadPool::get();
  const int nz = VElems.nz(), ny = VElems.ny();
  const int j_start = 1, j_end = ny - 1;
  const int k_start = 1, k_end = nz - 1;
  const int j_size = j_end - j_start;
  const int k_size = k_end - k_start;
  const size_t total_iterations = static_cast<size_t>(k_size * j_size);
  auto nChanges_futures = pool.submit_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) -> size_t {
        size_t local_nChanges = 0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          int j = j_start + (idx % j_size);
          int k = k_start + (idx / j_size);
          const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
          for (short ix = 0; ix < s.cnt; ++ix) {
            for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i) {
              // voxel* v=s.s[ix].v(i-1);
              // if (v)

              const int *pijk = &voxls(i, j, k);
              long pID = *pijk;

              if (pID == rawValue) {
                //			 float R = cg.segs_[(k - 1) * cg.ny + (j
                //- 1)].vxl(i - 1)->R;

                short nDifferentID = 0;
                if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk))
                  nDifferentID++;
                if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk))
                  nDifferentID++;
                if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk))
                  nDifferentID++;
                if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk))
                  nDifferentID++;
                if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk))
                  nDifferentID++;
                if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk))
                  nDifferentID++;

                if (nDifferentID >= 2) {
                  map<int, short> neis;

                  long neI = voxls.v_i(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_i(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_j(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_j(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_k(-1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);
                  neI = voxls.v_k(1, pijk);
                  if (neI != pID && bgn <= neI && lst >= neI)
                    ++(neis.insert(pair<int, short>(neI, 0)).first->second);

                  map<int, short>::iterator neitr =
                      max_element(neis.begin(), neis.end(), mapComparer<int>());
                  if (neitr->second >= 2) {
                    ++local_nChanges;
                    VElems(i, j, k) = neitr->first;
                  }
                }
              }
            }
          }
        }
        return local_nChanges;
      });
  size_t nChanges(0);
  for (auto &future : nChanges_futures) {
    nChanges += future.get();
  }

  cout << "  ngMedLoose:" << nChanges << "  ";
}

void medianElem(const inputDataNE &cg, voxelField<int> &VElems,
                voxelField<int> &voxls, long bgn, long lst,
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
          const size_t k = coords[0], j = coords[1];
          if (k == 0 || k == nz - 1 || j == 0 || j == ny - 1)
            continue;
          const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
          for (short ix = 0; ix < s.cnt; ++ix) {
            for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i) {
              const int *pijk = &voxls(i, j, k);
              const long pID = *pijk;

              if (pID < bgn || pID > lst)
                continue;

              short nSameID = 0;
              short nDifferentID = 0;
              std::array<int, 6> nei_buf{};
              short nei_count = 0;

              auto visit_neighbor = [&](long neID) {
                if (neID == pID) {
                  ++nSameID;
                } else if (neID >= bgn && neID <= lst) {
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

              if ((nDifferentID <= nSameID) || (nei_count == 0))
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

#endif
