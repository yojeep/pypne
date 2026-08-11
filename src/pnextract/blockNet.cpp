#include "blockNet.hpp"
#include "ElementGNE.hpp"
#include "globals.hpp"
#include "gtl/phmap.hpp"
#include <array>
#include <boost/sort/sort.hpp>
#include <cmath>
#include <iostream>
#include <vector>

void blockNetwork::CreateVElem() { /// ### map pore labels from maximal spheres
                                   /// to the image
  /// blockNetwork::VElems, standing for VoxelElements.
  BitMap3D<voxel> &vxlMap = srf.vxlMap;
  std::vector<medialBall> &ballSpace = srf.ballSpace;
  const BitVector3D &bool_img = srf.bool_img;

  poreAdditBalls.reserve(nBP6);

  std::cout << "\nCreating pore elements:" << std::endl;

  if (poreIs.empty()) {
    sides[0] = voxel(-1, cfg.ny / 2, cfg.nz / 2, _SideRadius); // minX
    poreAdditBalls.emplace_back(&sides[0], -1);
    poreIs.emplace_back(0, cfg.ny * cfg.nz, &poreAdditBalls.back());

    sides[1] = voxel(cfg.nx + 1, cfg.ny / 2, cfg.nz / 2, _SideRadius); // maxX
    poreAdditBalls.emplace_back(&sides[1], -1);
    poreIs.emplace_back(0, cfg.ny * cfg.nz, &poreAdditBalls.back());

    if (nBP6 == 6) {
      sides[2] = voxel(cfg.nx / 2, -1, cfg.nz / 2, _SideRadius); // minY
      poreAdditBalls.emplace_back(&sides[2], -1);
      poreIs.emplace_back(0, cfg.nx * cfg.nz, &poreAdditBalls.back());

      sides[3] = voxel(cfg.nx / 2, cfg.ny + 1, cfg.nz / 2, _SideRadius); // maxY
      poreAdditBalls.emplace_back(&sides[3], -1);
      poreIs.emplace_back(0, cfg.nx * cfg.nz, &poreAdditBalls.back());

      sides[4] = voxel(cfg.nx / 2, cfg.ny / 2, -1, _SideRadius); // minZ
      poreAdditBalls.emplace_back(&sides[4], -1);
      poreIs.emplace_back(0, cfg.ny * cfg.nx, &poreAdditBalls.back());

      sides[5] = voxel(cfg.nx / 2, cfg.ny / 2, cfg.nz + 1, _SideRadius); // maxZ
      poreAdditBalls.emplace_back(&sides[5], -1);
      poreIs.emplace_back(0, cfg.ny * cfg.nx, &poreAdditBalls.back());
    }
  }

  int unassigned = -1;
  VElems = Eigen::Tensor<int, 3, Eigen::RowMajor>(cfg.nz + 2, cfg.ny + 2,
                                                  cfg.nx + 2);
  VElems.setConstant(unassigned);
  VElems.chip(0, 2).setConstant(0);
  VElems.chip(cfg.nx + 1, 2).setConstant(1);
  int nVVs = 1;

  if (nBP6 == 2) {
    VElems.chip(0, 1).setConstant(-1 - 2 - nVVs);
    VElems.chip(cfg.ny + 1, 1).setConstant(-1 - 3 - nVVs);
    VElems.chip(0, 0).setConstant(-1 - 4 - nVVs);
    VElems.chip(cfg.nz + 1, 0).setConstant(-1 - 5 - nVVs);
  } else if (nBP6 == 6) {
    VElems.chip(0, 1).setConstant(2);
    VElems.chip(cfg.ny + 1, 1).setConstant(3);
    VElems.chip(0, 0).setConstant(4);
    VElems.chip(cfg.nz + 1, 0).setConstant(5);
  } else {
    std::cout << "nBP6 must be 2 or 6" << std::endl;
    exit(1);
  }
  // std::cout << "VElems:\n" << VElems << std::endl;
  std::cout << "Creating initial pore elements" << std::endl;

  {

    firstPore = (poreIs.size());
    for (const auto &bi : ballSpace) {
      if (bi.boss == &bi) {
        int iz = static_cast<int>(bi.fz) + 1;
        int iy = static_cast<int>(bi.fy) + 1;
        int ix = static_cast<int>(bi.fx) + 1;
        VElems(iz, iy, ix) = poreIs.size();
        poreIs.emplace_back();
        poreIs.back().mb = &bi;
      }
    }
    std::cout << "\n created " << poreIs.size()
              << " pores (+boundaries)  for up to index " << 0 << std::endl;
    lastPore = (poreIs.size() - 1);

    // const int firstPoreInd = firstPores;
    // const int lastPoreInd = lastPores;
    std::cout << " mapping pores indices to image, for index " << 0 << ":  "
              << firstPore << " to " << lastPore
              << ",  unassigned:" << unassigned << std::endl;

    for (const auto &bi : ballSpace) {
      if (!bi.boss) {
        std::cout << "skipping " << bi.fz << " " << bi.fy << " " << bi.fx
                  << std::endl;
        continue;
      }

      const medialBall *mastrSphere = bi.mastrSphere();
      int izm = static_cast<int>(mastrSphere->fz) + 1;
      int iym = static_cast<int>(mastrSphere->fy) + 1;
      int ixm = static_cast<int>(mastrSphere->fx) + 1;

      int Vm = VElems(izm, iym, ixm);

      float rbi = bi.R;
      float r = std::max(rbi * 0.5f - 1.0f, 1.001f);
      for_each_voxel_in_sphere(bi.fz, bi.fy, bi.fx, r, cfg.nz, cfg.ny, cfg.nx,
                               [&](int z, int y, int x) {
                                 if (!bool_img(z, y, x))
                                   return;
                                 int zVE = z + 1;
                                 int yVE = y + 1;
                                 int xVE = x + 1;
                                 int Vn = VElems(zVE, yVE, xVE);
                                 if (Vn == unassigned) {
                                   VElems(zVE, yVE, xVE) = Vm;
                                   return;
                                 }

                                 if (Vn == Vm)
                                   return;

                                 const voxel &vn = vxlMap(z, y, x);
                                 if (vn.ball || vn.R >= rbi)
                                   return;

                                 const medialBall *bn = poreIs[Vn].mb;
                                 float dold = sq(x + 0.5f - bn->fx) +
                                              sq(y + 0.5f - bn->fy) +
                                              sq(z + 0.5f - bn->fz);
                                 float dnew = sq(x + 0.5f - mastrSphere->fx) +
                                              sq(y + 0.5f - mastrSphere->fy) +
                                              sq(z + 0.5f - mastrSphere->fz);

                                 if (dnew >= dold)
                                   return;
                                 VElems(zVE, yVE, xVE) = Vm;
                               });
    }

    std::cout << "\n growing pores, ..  ";
    std::cout.flush();
    PoreGrower pg(VElems, vxlMap, ballSpace, firstPore, unassigned);
    pg.grow();
  }
}
// friend size_t hash_value(const Person &p)
// {
//     return phmap::HashState().combine(0, p._first, p._last, p._age);
// }
void blockNetwork::createNewThroats() {

  BitMap3D<voxel> &vxlMap = srf.vxlMap;

  struct hash_array {
    size_t operator()(const std::array<int, 2> &arr) const {
      return gtl::HashState().combine(0, arr[0], arr[1]);
    }
  };

  gtl::flat_hash_map<std::array<int, 2>, int, hash_array> unique_paires;

  {
    std::cout << "\nlooking for connections, iFirst = " << throatIs.size()
              << " ..  ";
    std::cout.flush();

    // Use fixed-size arrays to avoid any dynamic memory allocation
    const int dx[3] = {1, 0, 0};
    const int dy[3] = {0, 1, 0};
    const int dz[3] = {0, 0, 1};
    const int signs[6] = {1, 1, 1, -1, -1, -1};
    const int axes[6] = {0, 1, 2, 0, 1, 2};

    for (const voxel &vi : vxlMap.vec_data_) {
      int z0 = vi.z + 1;
      int y0 = vi.y + 1;
      int x0 = vi.x + 1;

      int p0ID = VElems(z0, y0, x0);
      if (p0ID == -1)
        continue;

      poreNE &p0 = poreIs[p0ID];
      if (p0ID >= firstPore)
        ++(p0.volumn);

      // Core trick: use a bitmask to combine directions to check
      // Bits 0,1,2 correspond to +X, +Y, +Z respectively.
      // By default only positive directions are checked (1|2|4 = 7).
      // If on a boundary (e.g. x0==1), also add the negative direction.
      int mask = 1 | 2 | 4;
      if (x0 == 1)
        mask |= 8; // add -X
      if (y0 == 1)
        mask |= 16; // add -Y
      if (z0 == 1)
        mask |= 32; // add -Z

      // Loop over all six possible directions (0..5)
      for (int d = 0; d < 6; ++d) {
        // Skip direction if not enabled by the mask
        if (!(mask & (1 << d)))
          continue;

        // Determine offset based on direction index
        int sign = signs[d];
        int axis = axes[d];

        int xn = x0 + dx[axis] * sign;
        int yn = y0 + dy[axis] * sign;
        int zn = z0 + dz[axis] * sign;

        int p1ID = VElems(zn, yn, xn);
        if (p1ID < 0 || p0ID == p1ID)
          continue;

        // Canonical ordering to ensure undirected graph connectivity is unique
        bool dir = p0ID < p1ID;
        int lID = dir ? p0ID : p1ID;
        int hID = dir ? p1ID : p0ID;

        poreNE &pLow = poreIs[lID];
        poreNE &pHigh = poreIs[hID];

        auto [it, inserted] = unique_paires.try_emplace(
            std::array<int, 2>{lID, hID}, static_cast<int>(throatIs.size()));

        if (inserted) {
          throatIs.emplace_back(0, lID, hID);
        }

        // Accumulate cross-sectional area contribution
        // Positive direction: +1, negative direction: -1
        throatNE &trot = throatIs[it->second];
        trot.CrosArea[axis] += dir ? sign : -sign;

        // Accumulate surface area for both pores
        if (lID >= firstPore)
          ++(pLow.surfaceArea);
        if (hID >= firstPore)
          ++(pHigh.surfaceArea);
      }
    }

    boost::sort::sample_sort(
        throatIs.begin(), throatIs.end(), [](const auto &a, const auto &b) {
          return a.e1 < b.e1 || (a.e1 == b.e1 && a.e2 < b.e2);
        });

    for (size_t tid = 0; tid < throatIs.size(); ++tid) {
      throatNE &elemt = throatIs[tid];
      elemt.tid = tid;
      poreNE &elemp1 = poreIs[elemt.e1];
      poreNE &elemp2 = poreIs[elemt.e2];
      elemp1.contacts.emplace(elemt.e2, tid);
      elemp2.contacts.emplace(elemt.e1, tid);
    }

    (std::cout << throatIs.size() << ", ").flush();
  }

  nNodes = poreIs.size();
  nTrots = throatIs.size();
  nElems = nNodes + nTrots;

  std::cout << "nElems:  " << nElems << " = " << poreIs.size() << " + "
            << throatIs.size() << std::endl;

  throadAdditBalls.reserve(
      throatIs.size() * 2); // to improve the efficiency when later generating
                            // and storing new maximal balls for throat surfaces

  std::cout << "\ncalc throat properties: ";
  std::cout.flush();

  {
    std::cout << " collecting face voxels,  ";
    std::cout.flush();

    for (auto &tr : throatIs) {
      size_t n_neis = std::abs(tr.CrosArea[0]) + std::abs(tr.CrosArea[1]) +
                      std::abs(tr.CrosArea[2]) + 1;
      tr.toxels1.reserve(n_neis);
      tr.toxels2.reserve(n_neis);
    }

    int nMultiTouchErrors = 0;

    for (const voxel &vi : vxlMap.vec_data_) {
      int z = vi.z, y = vi.y, x = vi.x;
      int zVE = z + 1;
      int yVE = y + 1;
      int xVE = x + 1;

      int p0ID = VElems(zVE, yVE, xVE);
      if (p0ID == -1)
        continue;
      gtl::flat_hash_set<int> neis;
      neis.reserve(6);

      auto add_nei = [&](int neiPID) noexcept {
        if ((p0ID != neiPID) && (neiPID >= 0))
          neis.emplace(neiPID);
      };
      add_nei(VElems(zVE, yVE, xVE - 1));
      add_nei(VElems(zVE, yVE, xVE + 1));
      add_nei(VElems(zVE, yVE - 1, xVE));
      add_nei(VElems(zVE, yVE + 1, xVE));
      add_nei(VElems(zVE - 1, yVE, xVE));
      add_nei(VElems(zVE + 1, yVE, xVE));
      poreNE &p0 = poreIs[p0ID];
      for (int nei : neis) {
        throatNE &trot = throatIs[p0.contacts[nei]];
        if (p0ID > nei) {
          trot.toxels2.emplace_back(&vxlMap(z, y, x));
        } else {
          trot.toxels1.emplace_back(&vxlMap(z, y, x));
        }
      }
      if (neis.size() > 1)
        ++nMultiTouchErrors;
    }

    if (nMultiTouchErrors > 0)
      std::cout << "\n   Warning: " << nMultiTouchErrors
                << " voxels being in touch to more than two pores" << std::endl;

    std::cout << " calculating throat radii";
    std::cout.flush();

    for (throatNE &tr : throatIs) {
      bool t1_empty = tr.toxels1.empty();
      bool t2_empty = tr.toxels2.empty();

      if (t1_empty && t2_empty) {
        std::cout << "  ERROR1017, toxls size: 0 0  ";
        std::cout.flush();
        continue;
      }

      if (!t1_empty)
        boost::sort::sample_sort(tr.toxels1.begin(), tr.toxels1.end(),
                                 metaballcomparer(), num_workers);
      if (!t2_empty)
        boost::sort::sample_sort(tr.toxels2.begin(), tr.toxels2.end(),
                                 metaballcomparer(), num_workers);

      auto processToxel = [&](std::vector<voxel *> &toxels, int type, int e_val,
                              int ball_param) noexcept {
        if (toxels.empty())
          return;

        voxel *tvox = toxels.front();
        medialBall *vbi = tvox->ball;

        if (vbi != nullptr) {
          vbi->type = type;
          medialBall *mvi = vbi->mastrSphere();
          if (mvi && mvi != vbi &&
              e_val != VElems(mvi->iz() + 1, mvi->iy() + 1, mvi->ix() + 1)) {
            std::cout << " Dmb" << (type == 6 ? "1" : "2") << "rrr  "
                      << VElems(mvi->iz() + 1, mvi->iy() + 1, mvi->ix() + 1)
                      << "   ";
          }
        } else {
          throadAdditBalls.emplace_back(tvox, ball_param);
          tvox->ball = &throadAdditBalls.back();
          srf.moveUphill(*tvox->ball);
        }
      };

      processToxel(tr.toxels1, 6, tr.e1, 16);
      processToxel(tr.toxels2, 5, tr.e2, 15);
    }

    std::cout << "." << std::endl;
  }
}