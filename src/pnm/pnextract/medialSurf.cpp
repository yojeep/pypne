
#include "medialSurf.h"
#include "ElementGNE.h"
#include "edt_double.h"
#include "globals.h"
#include "inputData.h"
// #include "typses.h"
#include "indexunraveler.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <boost/sort/sort.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>

medialSurface::medialSurface(
    inputDataNE &cfg) //, double vmvLimRelF, double crossAreaf
    : cg_(cfg), segs_(cfg.segs_),
      binary_image(cfg.VImage.data_.data(), cfg.nz, cfg.ny, cfg.nx),
      binary_image_flat(cfg.VImage.data_.data(), cfg.nz * cfg.ny * cfg.nx) {
  setDefaults(
      -5.); // set _minRp to negative value if not provided by user, to be
            // re-assigned in calc_distmaps(),  to be synced with setDefaults()

  nx = cfg.nx;
  ny = cfg.ny;
  nz = cfg.nz;
  auto &pool = GlobalThreadPool::get();
  const int total_iterations = nz * ny;
  std::atomic<size_t> _nVxls(0);
  pool.detach_blocks(0, total_iterations,
                     [&](const size_t start, const size_t end) {
                       size_t local_nVxls = 0;
                       for (size_t i = start; i < end; ++i) {
                         const segments &s = cg_.segs_[i];
                         for (int ix = 0; ix < s.cnt; ++ix)
                           if (s.s[ix].value == 0)
                             local_nVxls += s.s[ix + 1].start - s.s[ix].start;
                       }
                       _nVxls.fetch_add(local_nVxls, std::memory_order_relaxed);
                     });
  pool.wait();
  nVxls = _nVxls.load();
  invalidSeg.start = -10000;
  invalidSeg.value = 255;
}

void medialSurface::setDefaults(double avgR) {
  /// minRPore/Rnoise is a different keyword, of its own, not part of
  /// medialSurfaceSettings.    It is supposed to be the only adjustable
  /// parameter for the regular users.       The default value of minRPore
  /// is 1.75.

  /// medialSurfaceSettings are for advanced users and its defaults are chosen
  /// (if not provided by the user) based on the minRPore value.

  /// To extract a network that has a lower network coordination number, you can
  /// decrease the values of lenNf (e.g. 0.4), and vmvRadRelNf (e.g. 1.05) and
  /// increase the values of nRSmoothing (e.g. 9), RCorsf (e.g. 0.2) and RCors
  /// (e.g. 2.5).  You can also consider increasing minRPore (also named Rnoise)
  /// keyword to let say 2..     Every change you make you need to check that
  /// the network produces reasonable results as these are sensitive parameters
  /// and do not behave linearly. We should not be woried about the high
  /// coordination number as long as the network predicts the physical
  /// properties correctly, but not everybody agrees with me here!

  _minRp = std::min(1.25, avgR * 0.25) + 0.5;
  if (cg_.giv("Rnoise" + _s(0), _minRp) || cg_.giv("minRPore", _minRp) ||
      cg_.giv("Rnoise", _minRp))
    std::cout << " minimum pore radius: " << _minRp << std::endl;
  else
    std::cout << " keyword \"minRPore\" not found, default value ("
              << std::abs(_minRp) << ") will be used" << std::endl;

  _clipROutx = 0.05;
  _clipROutyz = 0.98;
  _midRf = 0.7;
  _MSNoise = 1. * std::abs(_minRp) + 1.;
  _lenNf = 0.6;
  _vmvRadRelNf = 1.1;
  _nRSmoothing = 3;
  _RCorsnf = 0.15;
  _RCorsn = std::abs(_minRp);

  if (cg_.nBP6 == 6)
    _clipROutyz = _clipROutx;

  std::istringstream keywrdData;
  if (cg_.giv("medialSurfaceSettings" + _s(0), keywrdData) ||
      cg_.giv("medialSurfaceSettings", keywrdData)) {
    keywrdData >> _clipROutx >> _clipROutyz >> _midRf >> _MSNoise >> _lenNf >>
        _vmvRadRelNf >> _nRSmoothing >> _RCorsnf >> _RCorsn;
  }

  if (_minRp < 0.)
    std::cout << " Default setting, will be updated after distance map "
                 "computation:\n";
  std::cout << "  minRPore     : " << abs(_minRp) << ";\n";
  std::cout << "  medialSurfaceSettings: " << _clipROutx << "  " << _clipROutyz
            << "  " << _midRf << "  " << _MSNoise << "  " << _lenNf << "  "
            << _vmvRadRelNf << "  " << _nRSmoothing << "  " << _RCorsnf << "  "
            << _RCorsn << std::endl;

  std::cout << "  medialSurfaceSettings:\n"
            << "   clipROutx     : " << _clipROutx << "\n"
            << "   clipROutyz    : " << _clipROutyz << "\n"
            << "   midRFrac      : " << _midRf << "\n"
            << "   RMedSurfNoise : " << _MSNoise << "\n"
            << "   lenNf         : " << _lenNf << "\n"
            << "   vmvRadRelNf   : " << _vmvRadRelNf << "\n"
            << "   nRSmoothing   : " << _nRSmoothing << "\n"
            << "   RCorsf  : " << _RCorsnf << "\n"
            << "   RCors   : " << _RCorsn << "\n"
            << std::endl;
}

template <typename T> union ndi_t {
  T value;

  // 关键：默认构造函数什么都不做
  ndi_t() {}

  // 允许显式构造
  template <typename... Args>
  ndi_t(Args &&...args) : value(static_cast<Args &&>(args)...) {}

  // 隐式转换成 T&
  operator T &() { return value; }
  T *operator->() { return &value; }
};

void medialSurface::buildvoxelspace() { ///  Build voxelspace -- memory for
                                        ///  void/active voxels
  std::cout << "\nProcessing " << cg_._rockTypes[0].name
            << " voxels:" << std::endl;
  std::cout << " Creating " << nVxls << " voxels with index: " << int(0);
  std::cout.flush();
  // binary_image_flat =
  //     std::span<const uint8_t>(cg_.VImage.data_.data(), nz * ny * nx);

  // binary_image = std::mdspan<const uint8_t, std::dextents<size_t, 3>>(
  //     cg_.VImage.data_.data(), nz, ny, nx);

  // vxlSpace = Eigen::array<voxel, 1>(int(nVxls));

  vxlSpace.reserve(nVxls);
  vxlMap = Eigen::Tensor<voxel *, 3, Eigen::RowMajor>(nz, ny, nx);

  IndexUnraveler unraveler({
      static_cast<size_t>(nz),
      static_cast<size_t>(ny),
      static_cast<size_t>(nx),
  });
  size_t total_iterations = nz * ny * nx;

  auto coords = std::array<int, 3>{};
  size_t count = 0;
  for (size_t idx = 0; idx < total_iterations; ++idx) {
    unraveler.unravel(idx, coords);
    int z = coords[0], y = coords[1], x = coords[2];
    if (binary_image_flat(idx) == 0) {
      vxlSpace.emplace_back(x, y, z, 0);
      voxel &vi = vxlSpace[count];
      vi.i = x;
      vi.j = y;
      vi.k = z;
      vxlMap(z, y, x) = &vi;
      ++count;
    } else {
      vxlMap(z, y, x) = nullptr;
    }
  }

  // for (int iz = 0; iz < nz; ++iz) {
  //   for (int iy = 0; iy < ny; ++iy) {
  //     const segments &s = segs_[iz * ny + iy];
  //     for (int ix = 0; ix < s.cnt; ++ix) {
  //       if (s.s[ix].value == 0)
  //         for (int i = s.s[ix].start; i < s.s[ix + 1].start; ++i) {
  //           p->i = i;
  //           p->j = iy;
  //           p->k = iz;
  //           // size_t currentIndex = p - vxlBegin;
  //           // iZ[iz].push_back(currentIndex);
  //           ++p;
  //         }
  //     }
  //   }
  // }

  /// Link voxels to segments
  auto p = vxlSpace.begin();
  total_iterations = nz * ny;
  for (size_t izy = 0; izy < total_iterations; ++izy) {
    segments &s = segs_[izy];
    for (int ix = 0; ix < s.cnt; ++ix)
      if (s.s[ix].value == 0) {
        s.s[ix].segV = &*p;
        p += s.s[ix + 1].start - s.s[ix].start;
      }
  }
}

void medialSurface::paradox_pre_removeincludedballI() // to remove the included
                                                      // maximal bals
{ /// Remove maximal-balls, leave one in each adjacent voxel. This saves time
  /// when sorting in paradoxremoveincludedballI()
  if (!nVxls) {
    return;
  }
  std::cout << " pre-remove included balls: out of " << vxlSpace.size();
  std::cout.flush();

  auto &pool = GlobalThreadPool::get();
  const size_t z_steps = (nz + 1) / 2;
  const size_t y_steps = (ny + 1) / 2;
  const size_t x_steps = (nx + 1) / 2;
  const size_t total_iterations = z_steps * y_steps * x_steps;

  IndexUnraveler unraveler({z_steps, y_steps, x_steps});

  std::atomic<uint64_t> _nBalls(0);
  pool.detach_blocks(
      0, total_iterations, [&](const size_t start, const size_t end) {
        std::array<int, 3> coord{};
        for (size_t idx = start; idx < end; ++idx) {
          unraveler.unravel(idx, coord);
          int zo = coord[0] * 2, yo = coord[1] * 2, xo = coord[2] * 2;
          float max_r = -1.0;
          voxel *max_voxel = nullptr;
          for (short zi = 0; zi < 2; ++zi)
            for (short yi = 0; yi < 2; ++yi)
              for (short xi = 0; xi < 2; ++xi) {
                int z = zo + zi, y = yo + yi, x = xo + xi;
                if (z < 0 || z >= nz || y < 0 || y >= ny || x < 0 || x >= nx ||
                    binary_image(z, y, x) != 0)
                  continue;
                voxel *v = vxlMap(z, y, x);

                float v_r = v->R;
                if (v_r > max_r && v_r > _minRp) {
                  max_r = v_r;
                  max_voxel = v;
                }
              }
          if (max_voxel != nullptr) {
            max_voxel->ball = &ToBeAssigned;
            _nBalls.fetch_add(1, std::memory_order_relaxed);
          }
        }
      });

  pool.wait();
  nBalls = _nBalls.load();
}

void medialSurface::paradoxremoveincludedballI() { /// Remove included balls.
                                                   /// What remains are called
                                                   /// maximal spheres
  if (!nVxls) {
    return;
  }

  std::cout << " sorting... ";
  std::vector<voxel *> tvs;
  tvs.reserve(nBalls);

  for (voxel &v : vxlSpace) {
    if (v.ball) {
      tvs.push_back(&v);
    }
  }

  std::cout << tvs.size() << " balls" << std::endl;
  boost::sort::sample_sort(
      tvs.begin(), tvs.end(),
      [](const voxel *a, const voxel *b) { return a->R > b->R; }, num_workers);

  std::cout << " remove included balls:";
  std::cout.flush();

  size_t ndel = 0;
  for (size_t idx = 0; idx < tvs.size(); ++idx) {
    voxel *vi = tvs[idx];
    if (!vi->ball)
      continue;

    int zo = vi->k, yo = vi->j, xo = vi->i;
    float ri = vi->R;
    float ripinc = ri + 0.55; //.+RPreDelete
    float mbmbDist = _RCorsnf * ri + _RCorsn;
    float ripinc2 = ripinc * ripinc;
    int rz = ripinc;

    for (int zi = -rz; zi < rz + 1; ++zi) {
      float ry2 = ripinc2 - zi * zi;
      if (ry2 <= 0)
        continue;
      int ry = std::sqrtf(ry2);
      for (int yi = -ry; yi < ry + 1; ++yi) {
        float rx2 = ry2 - yi * yi;
        if (rx2 <= 0)
          continue;
        int rx = std::sqrtf(rx2); // sqrts(r2i)+1-a-b;
        for (int xi = -rx; xi < rx + 1; ++xi) {
          int z = zo + zi, y = yo + yi, x = xo + xi;
          if (z < 0 || z >= nz || y < 0 || y >= ny || x < 0 || x >= nx ||
              binary_image(z, y, x) != 0)
            continue;
          voxel *vj = vxlMap(z, y, x);
          if ((vj->ball) && (vi != vj)) {
            float rj = vj->R;
            if (rj <= ri) {
              float D = std::sqrtf(zi * zi + yi * yi + xi * xi);
              if (D < mbmbDist || (D + rj < ripinc + _MSNoise)) {
                vj->ball = nullptr;
                ++ndel;
              }
            }
          }
        }
      }
    }

    if ((idx % 10000) == 0)
      std::cout << "\r  remove = " << ndel;
  }
  std::cout << "\r  removed = " << ndel << " remained = " << tvs.size() - ndel
            << " balls" << std::endl;
  nBalls -= ndel;
  return;
}

void medialSurface::moveUphill(medialBall *b_i) // const
{ /// Refines the maximal-sphere location and radius,

  const voxel *vi = vxlMap(b_i->iz(), b_i->iy(), b_i->ix());
  const int zm = vi->k, ym = vi->j, xm = vi->i;
  const int zn = zm - 1, zp = zm + 1, yn = ym - 1, yp = ym + 1, xn = xm - 1,
            xp = xm + 1;
  std::array<float, 3> disp{};
  if (0 <= xn && xp < nx) {
    const voxel *vjn = vxlMap(zm, ym, xn);
    const voxel *vjp = vxlMap(zm, ym, xp);
    if (vjn && vjp) {
      float gp = vjp->R - vi->R;
      float gm = vi->R - vjn->R;
      if (std::abs(gp - gm) > 0.01)
        disp[0] = std::max(-0.49, std::min(0.49, -0.5 * (gp + gm) / (gp - gm)));
    }
  }
  if (0 <= yn && yp < ny) {
    const voxel *vjn = vxlMap(zm, yn, xm);
    const voxel *vjp = vxlMap(zm, yp, xm);
    if (vjn && vjp) {
      float gp = vjp->R - vi->R;
      float gm = vi->R - vjn->R;
      if (std::abs(gp - gm) > 0.01)
        disp[1] = std::max(-0.49, std::min(0.49, -0.5 * (gp + gm) / (gp - gm)));
    }
  }
  if (0 <= zn && zp < nz) {
    const voxel *vjn = vxlMap(zn, ym, xm);
    const voxel *vjp = vxlMap(zp, ym, xm);
    if (vjn && vjp) {
      float gp = vjp->R - vi->R;
      float gm = vi->R - vjn->R;
      if (std::abs(gp - gm) > 0.01)
        disp[2] = std::max(-0.49, std::min(0.49, -0.5 * (gp + gm) / (gp - gm)));
    }
  }
  if (b_i != b_i->boss) {
    std::array<float, 3> BosKidVec = {
        b_i->fi - b_i->boss->fi, // x
        b_i->fj - b_i->boss->fj, // y
        b_i->fk - b_i->boss->fk  // z
    };

    float dot = BosKidVec[0] * disp[0] + BosKidVec[1] * disp[1] +
                BosKidVec[2] * disp[2];

    // magSqr(BosKidVec)
    float magSq = BosKidVec[0] * BosKidVec[0] + BosKidVec[1] * BosKidVec[1] +
                  BosKidVec[2] * BosKidVec[2];

    // disp -= 0.95 * (dot / (magSq + 1e-12)) * BosKidVec
    float coeff = 0.95 * dot / (magSq + 1e-12);
    disp[0] -= coeff * BosKidVec[0];
    disp[1] -= coeff * BosKidVec[1];
    disp[2] -= coeff * BosKidVec[2];
  }
  b_i->fi = vi->i - _mp5 + disp[0];
  b_i->fj = vi->j - _mp5 + disp[1];
  b_i->fk = vi->k - _mp5 + disp[2];
  b_i->R = vi->R + 0.95 * std::sqrtf(disp[0] * disp[0] + disp[1] * disp[1] +
                                     disp[2] * disp[2]);
}

void medialSurface::moveUphillp1(medialBall *bi) // const
{ /// Refines the maximal-sphere location, by moving it uphil the gradient of
  /// the distance-map, potentially relocates to new voxels

  const voxel *vi = vxlMap(bi->iz(), bi->iy(), bi->ix());
  const int zm = vi->k, ym = vi->j, xm = vi->i;
  const int zn = zm - 1, zp = zm + 1, yn = ym - 1, yp = ym + 1, xn = xm - 1,
            xp = xm + 1;
  std::array<float, 3> disp{}, grad{};

  if (0 <= xn && xp < nx) {
    const voxel *vjm = vxlMap(zm, ym, xn);
    const voxel *vjp = vxlMap(zm, ym, xp);
    float gp = vjp->R - vi->R;
    float gm = vi->R - vjm->R;
    grad[0] = 0.5 * (gp + gm);
    if (std::abs(gp - gm) > 0.01)
      disp[0] = std::max(-0.59, std::min(0.59, -0.5 * (gp + gm) / (gp - gm)));
  }

  if (0 <= yn && yp < ny) {
    const voxel *vjm = vxlMap(zm, yn, xm);
    const voxel *vjp = vxlMap(zm, yp, xm);
    float gp = vjp->R - vi->R;
    float gm = vi->R - vjm->R;
    grad[1] = 0.5 * (gp + gm);
    if (std::abs(gp - gm) > 0.01)
      disp[1] = std::max(-0.59, std::min(0.59, -0.5 * (gp + gm) / (gp - gm)));
  }

  if (0 <= zn && zp < nz) {
    const voxel *vjm = vxlMap(zn, ym, xm);
    const voxel *vjp = vxlMap(zp, ym, xm);
    float gp = vjp->R - vi->R;
    float gm = vi->R - vjm->R;
    grad[2] = 0.5 * (gp + gm);
    if (std::abs(gp - gm) > 0.01)
      disp[2] = std::max(-0.59, std::min(0.59, -0.5 * (gp + gm) / (gp - gm)));
  }

  for (int d = 0; d < 3; ++d)
    disp[d] += 1.4f * grad[d];

  if (bi != bi->boss) {
    std::array<float, 3> BosKidVec = {
        bi->fi - bi->boss->fi, // x
        bi->fj - bi->boss->fj, // y
        bi->fk - bi->boss->fk  // z
    };

    float dot = BosKidVec[0] * disp[0] + BosKidVec[1] * disp[1] +
                BosKidVec[2] * disp[2];

    // magSqr(BosKidVec)
    float magSq = BosKidVec[0] * BosKidVec[0] + BosKidVec[1] * BosKidVec[1] +
                  BosKidVec[2] * BosKidVec[2];

    // disp -= 0.95 * (dot / (magSq + 1e-12)) * BosKidVec
    float coeff = 0.95 * dot / (magSq + 1e-12);
    disp[0] -= coeff * BosKidVec[0];
    disp[1] -= coeff * BosKidVec[1];
    disp[2] -= coeff * BosKidVec[2];
  }

  float len = 0.0f;
  for (int d = 0; d < 3; ++d)
    len += disp[d] * disp[d];
  float scale = 1.0f / (0.55f * std::sqrtf(len) + 0.05f);
  for (int d = 0; d < 3; ++d)
    disp[d] *= scale;

  int iz_vj = bi->fk + disp[2], iy_vj = bi->fj + disp[1],
      ix_vj = bi->fi + disp[0];
  if (0 <= iz_vj && iz_vj < nz && 0 <= iy_vj && iy_vj < ny && 0 <= ix_vj &&
      ix_vj < nx) {
    voxel *vj = vxlMap(iz_vj, iy_vj, ix_vj);
    if (vj && vi != vj && vj->R > vi->R && vj->ball == nullptr) {
      bi->vxl->ball = nullptr;
      bi->vxl = vj;
      bi->fi = vj->i - _mp5;
      bi->fj = vj->j - _mp5;
      bi->fk = vj->k - _mp5;
      bi->R = vj->R;
      vj->ball = bi;
    }
  }
}

// void makeFriend(medialBall *vi, medialBall *vj) {
//   if (vj->R > vi->R) {
//     medialBall *tmp = vi;
//     vi = vj;
//     vj = tmp;
//   }
//   if ((vi->R < 1.5 * vj->R) && (!vi->isNei(vj)) && (!vi->inParents(vj)) &&
//       !vj->inParents(vi)) {
//     vi->addNei(vj);
//     vj->addNei(vi);
//   }
// }

/*inline double cosAngleWithBossPerD(const medialBall* a, const medialBall* b)
{ dbl3 v1=*a-*b; if (b==b->boss) return 1./(mag(v1)); dbl3 v2=*b-*(b->boss);
        double dotProd=v2&v1;
        return sqrt(dotProd*dotProd/(magSqr(v1)*magSqr(v1)*magSqr(v2)));
}*/

void medialSurface::competeForParent(medialBall *vi, medialBall *vj) {
  const double noise = _MSNoise;

  const double ri = vi->R;
  const double rj = vj->R;
  const double riSqr = ri * ri;
  const double rjSqr = rj * rj;
  const double dSqr = distSqr(vi, vj);
  const double dVal = sqrt(dSqr);

  const double wsinv = 1.0 / (riSqr + rjSqr);
  int iz_middle = wsinv * (vi->fk * rjSqr + vj->fk * riSqr),
      iy_middle = wsinv * (vi->fj * rjSqr + vj->fj * riSqr),
      ix_middel = wsinv * (vi->fi * rjSqr + vj->fi * riSqr);

  if (iz_middle < 0 || iz_middle >= nz || iy_middle < 0 || iy_middle >= ny ||
      ix_middel < 0 || ix_middel >= nx)
    return;

  const voxel *middlevxl = vxlMap(iz_middle, iy_middle, ix_middel);

  if (!middlevxl)
    return;

  const double minR = std::min(ri, rj);
  if (middlevxl->R <= minR * _midRf - 0.5)
    return;
  if (1.01 * dVal >= ri + rj + 1.0 + noise)
    return;

  // Handle boss assignment cases
  if (vj->boss == vj && vi->mastrSphere() != vj) {
    if (ri >= rj)
      vj->boss = vi;
    else if (vi->boss->R <= rj)
      vi->boss = vj;
    else if (ri >= rj - noise && ri * _vmvRadRelNf + noise >= rj)
      vj->boss = vi;
  } else if (vi->boss == vi && vj->mastrSphere() != vi) {
    if (rj >= ri)
      vi->boss = vj;
    else if (vj->boss->R <= ri)
      vj->boss = vi;
    else if (rj >= ri - noise && rj * _vmvRadRelNf + noise >= ri)
      vi->boss = vj;
  }

  medialBall *mvi = vi->mastrSphere();
  medialBall *mvj = vj->mastrSphere();

  if (mvi == vj || mvj == vi)
    return;

  if (mvi == mvj) {
    const short leveli = vi->level();
    const short levelj = vj->level();
    const short levelDiff = leveli - levelj;

    const double distViBoss = dist(vi->boss, vi);
    const double distVjBoss = dist(vj->boss, vj);
    const double distViVj = dist(vi, vj);

    if (levelDiff < -1) // leveli + 1 < levelj
    {
      if ((vj->boss->R - vj->R + 2 * noise) / (distVjBoss + 0.25) <
          (vi->R - vj->R + 2 * noise + 0.01) / (distViVj + 0.2))
        vj->boss = vi;
    } else if (levelDiff > 1) // leveli > levelj + 1
    {
      if ((vi->boss->R - vi->R + 2 * noise) / (distViBoss + 0.25) <
          (vj->R - vi->R + 2 * noise + 0.01) / (distViVj + 0.2))
        vi->boss = vj;
    } else {
      if (levelDiff > 0) // leveli > levelj
      {
        if ((vi->boss->R - vi->R + 2 * noise) / (distViBoss + 1.2) <
                (vj->R - vi->R + 2 * noise) / (distViVj + 1.3) &&
            !vj->inParents(vi))
          vi->boss = vj;
      } else if (levelDiff < 0) // leveli < levelj
      {
        if ((vj->boss->R - vj->R + 2 * noise) / (distVjBoss + 1.2) <
                (vi->R - vj->R + 2 * noise) / (distViVj + 1.3) &&
            !vi->inParents(vj))
          vj->boss = vi;
      }

      // if (middlevxl->R >= 0.45 * (ri + rj) - 1.0 && dVal < (ri + rj) * 0.5
      // + 2.0) 	makeFriend(vi, vj);
    }
  } else // mvi != mvj
  {
    const double avgR = 0.5 * (mvi->R + mvj->R);
    if (distSqr(mvi, mvj) <= _lenNf * (avgR + 2 * noise) * (avgR + 2 * noise)) {
      // Ensure mvi is the larger one
      if (mvi->R < mvj->R) {
        std::swap(vi, vj);
        std::swap(mvi, mvj);
      }

      if (mvj->R < _vmvRadRelNf * vj->R + noise &&
          mvj->R < _vmvRadRelNf * vi->R + noise &&
          mvj->R < _vmvRadRelNf * vi->boss->R + noise) {
        while (vj != vj->boss && mvj->R < _vmvRadRelNf * vj->boss->R + noise) {
          medialBall *pvj = vj->boss;
          vj->boss = vi;
          vi = vj;
          vj = pvj;
        }
        if (vj->boss == vj && vi->mastrSphere() != vj)
          vj->boss = vi;
      }
    }

    if (vi != vj->boss) {
      mvi = vi->mastrSphere();
      mvj = vj->mastrSphere();

      short leveli = vi->level();
      short levelj = vj->level();
      const double distMviMvj = dist(mvi, mvj);
      const double distAvg = distMviMvj + 0.5 * noise;

      while (leveli >= levelj) {
        const double viBossRatio =
            (vi->boss->R - vi->R + 0.55 * noise) / (dist(mvi, vi) + distAvg);
        const double vjRatio =
            (vj->R - vi->R + 0.5 * noise) / (dist(mvj, vi) + distAvg);
        if (viBossRatio >= vjRatio)
          break;

        medialBall *pvi = vi->boss;
        vi->boss = vj;
        vj = vi;
        vi = pvi;
        ++levelj;
        --leveli;
      }

      while (levelj >= leveli) {
        const double vjBossRatio =
            (vj->boss->R - vj->R + 0.55 * noise) / (dist(mvj, vj) + distAvg);
        const double viRatio =
            (vi->R - vj->R + 0.5 * noise) / (dist(mvi, vj) + distAvg);
        if (vjBossRatio >= viRatio)
          break;

        medialBall *pvj = vj->boss;
        vj->boss = vi;
        vi = vj;
        vj = pvj;
        ++leveli;
        --levelj;
      }

      // makeFriend(vi, vj);
    }
  }
}

void medialSurface::findBoss(medialBall *vi) {

  float zo = vi->fk, yo = vi->fj, xo = vi->fi;
  float ripp = vi->R * 0.6 + 2. * _MSNoise + 2.;
  float ripp2 = ripp * ripp;

  float rz = ripp;
  for (float zi = -rz; zi < rz + 1e-6f; zi += 1.f) {
    float ry2 = ripp2 - zi * zi;
    if (ry2 <= 0)
      continue;
    float ry = std::sqrtf(ry2);
    for (float yi = -ry; yi < ry + 1e-6f; yi += 1.0f) {
      float rx2 = ry2 - yi * yi;
      if (rx2 <= 0)
        continue;
      float rx = std::sqrtf(rx2);
      for (float xi = -rx; xi < rx + 1e-6f; xi += 1.0f) {
        int z = zo + zi, y = yo + yi, x = xo + xi;
        if (z < 0 || z >= nz || y < 0 || y >= ny || x < 0 || x >= nx ||
            binary_image(z, y, x) != 0)
          continue;
        voxel *vj = vxlMap(z, y, x);

        if (vj->ball && (vi != vj->ball)) {
          competeForParent(vi, vj->ball);
        }
      }
    }
  }
}

voxelImage segToVxlMesh(
    const medialSurface &ref) { /// converts segments back to voxelImage
  int nx = ref.nx, ny = ref.ny, nz = ref.nz;
  voxelImage vxls(nx, ny, nz, 255);
  auto &pool = GlobalThreadPool::get();
  const int total_iterations = ref.nz * ref.ny; // collapse(2)

  // 使用 detach_blocks 并行执行
  pool.detach_blocks(
      0, total_iterations, [&](const size_t start, const size_t end) {
        for (size_t iter = start; iter < end; ++iter) {
          int iz = iter / ny;
          int iy = iter % ny;
          const segments &s = ref.segs_[iter];
          for (int ix = 0; ix < s.cnt; ++ix) {
            std::fill(&vxls(s.s[ix].start, iy, iz),
                      &vxls(s.s[ix + 1].start, iy, iz), s.s[ix].value);
          }
        }
      });
  pool.wait();
  return vxls;
}

void medialSurface::calc_distmaps() // search  MBs at each voxel
{

  std::cout << " computing distance map for index " << int(0);
  std::cout.flush();

  if (!nVxls) {
    std::cout << " no voxels no balls,\n" << std::endl;
    exit(1);
  }

  auto &pool = GlobalThreadPool::get();
  auto &binary_image = cg_.VImage.data_;
  size_t nvoxels = binary_image.size();

  auto src = binary_image.data();
  auto inv_labels = std::make_unique_for_overwrite<bool[]>(nvoxels);
  pool.detach_blocks(0, nvoxels,
                     [&](const size_t start_idx, const size_t end_idx) {
                       for (size_t idx = start_idx; idx < end_idx; ++idx)
                         inv_labels[idx] = (src[idx] == 0);
                     });
  pool.wait();

  auto dt2 = std::make_unique_for_overwrite<float[]>(nvoxels);

  edt::binary_edtsq<bool>(inv_labels.get(), static_cast<int64_t>(cg_.nx),
                          static_cast<int64_t>(cg_.ny),
                          static_cast<int64_t>(cg_.nz), 1.0f, 1.0f, 1.0f, false,
                          num_workers, dt2.get());

  const double clipROutyz = _clipROutyz;
  const double clipROutx = _clipROutx;

  std::atomic<double> totalR(0.0);

  pool.detach_blocks(
      0, nVxls, [&](const size_t start_idx, const size_t end_idx) {
        double local_R = 0.0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          voxel &vi = vxlSpace[idx];
          const int z = vi.k, y = vi.j, x = vi.i;

          double limit =
              static_cast<double>(std::sqrt(dt2[z * ny * nx + y * nx + x])) -
              0.5;

          double iSqr = std::min(static_cast<double>(y + 2),
                                 static_cast<double>(ny - y + 1));
          if (iSqr < limit)
            limit =
                std::max((1.0 - clipROutyz) * limit + clipROutyz * iSqr, 0.01);

          iSqr = std::min(static_cast<double>(z + 2),
                          static_cast<double>(nz - z + 1));
          if (iSqr < limit)
            limit =
                std::max((1.0 - clipROutyz) * limit + clipROutyz * iSqr, 0.01);

          iSqr = std::min(static_cast<double>(x + 2),
                          static_cast<double>(nx - x + 1));
          if (iSqr < limit)
            limit = std::max((1.0 - clipROutx) * limit + clipROutx * iSqr, 0.1);

          vi.R = limit;
          local_R += limit;
        }
        totalR.fetch_add(local_R, std::memory_order_relaxed);
      });

  pool.wait();
  double rBalls = totalR.load();
  double avgR = rBalls / nVxls;

  std::cout << "\n";
  std::cout << "  average distance map = " << avgR << std::endl;

  if (_minRp < 0.) {
    setDefaults(avgR);
  }
}

void medialSurface::smoothRadius() {

  (std::cout << " smoothing R  ").flush();

  // std::vector<float> delRrr(vxlSpace.size(), 0.0f);
  (std::cout << "*").flush();

  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = nVxls;
  size_t nvoxels = binary_image_flat.size();
  if (!nvoxels) {
    std::cout << " no voxels no balls,\n" << std::endl;
    exit(1);
  }
  // auto _delRrr = std::make_unique_for_overwrite<float[]>(nvoxels);
  // std::mdspan<float, std::dextents<size_t, 3>> delRrr(_delRrr.get(), cg_.nz,
  //                                                     cg_.ny, cg_.nx);

  Eigen::Tensor<float, 3, Eigen::RowMajor> delRrr(cg_.nz, cg_.ny, cg_.nx);

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          voxel &vi = vxlSpace[idx];
          const int zo = vi.k, yo = vi.j, xo = vi.i;
          double sumR{0.};
          size_t counter{0};
          for (short zi = -1; zi < 2; ++zi)
            for (short yi = -1; yi < 2; ++yi)
              for (short xi = -1; xi < 2; ++xi) {
                int z = zo + zi, y = yo + yi, x = xo + xi;
                if (z < 0 || z >= nz || y < 0 || y >= ny || x < 0 || x >= nx ||
                    binary_image(z, y, x) != 0)
                  continue;
                sumR += vxlMap(z, y, x)->R;
                ++counter;
              }
          delRrr(zo, yo, xo) = 4. * sumR / (3 * counter + 27) - vi.R;
        }
      });
  pool.wait();
  (std::cout << "*").flush();

  pool.detach_blocks(
      0, total_iterations, [&](const size_t start_idx, const size_t end_idx) {
        for (size_t idx = start_idx; idx < end_idx; ++idx) {

          voxel &vi = vxlSpace[idx];
          const int zo = vi.k, yo = vi.j, xo = vi.i;
          double sumDelR{0.};
          size_t counter{0};
          for (short zi = -1; zi < 2; ++zi)
            for (short yi = -1; yi < 2; ++yi)
              for (short xi = -1; xi < 2; ++xi) {
                int z = zo + zi, y = yo + yi, x = xo + xi;
                if (z < 0 || z >= nz || y < 0 || y >= ny || x < 0 || x >= nx ||
                    binary_image(z, y, x) != 0)
                  continue;
                sumDelR += delRrr(z, y, x);
                ++counter;
              }
          vi.R +=
              std::min(std::max(0.02 * (delRrr(zo, yo, xo) -
                                        0.99 * 2. * sumDelR / (counter + 27)),
                                -0.005),
                       0.01);
        }
      });
  pool.wait();
  (std::cout << "*").flush();

  /// Finally, report max distance map, to confirm that distance map is not
  /// changed too much
  // float maxrrr = 0.0f;
  // for (size_t i = 0; i < vxlSpace.size(); ++i)
  //   maxrrr = std::max(maxrrr, vxlSpace[i].R);

  float maxrrr =
      std::max_element(vxlSpace.begin(), vxlSpace.end(),
                       [](const auto &a, const auto &b) { return a.R < b.R; })
          ->R;

  std::cout << " maxrrr " << maxrrr << std::endl;
}

void medialSurface::createBallsAndHierarchy() { /// Create distance map,
                                                /// maximal-spheres, and their
                                                /// hirarchy (medial-surface
                                                /// connectivity)
  //	mediaAxes medAxis(*this, minRP, clipOutSideBallFraction,
  // clipOutSideBallFraction*0.5+0.);

  buildvoxelspace();

  calc_distmaps();

  for (int i = 0; i < _nRSmoothing; ++i)
    smoothRadius();

  // std::atomic<size_t> _nBalls(0);
  // std::atomic<double> rBalls(0.0);

  auto &pool = GlobalThreadPool::get();

  paradox_pre_removeincludedballI();
  std::cout << "\n  number of potential maximal spheres: " << nBalls
            << std::endl;

  paradoxremoveincludedballI();

  std::cout << " collecting maximal balls out of " << nBalls << std::endl;

  std::vector<voxel *> tvs;
  tvs.reserve(nBalls);

  for (voxel &v : vxlSpace) {
    if (v.ball) {
      tvs.push_back(&v);
    }
  }

  boost::sort::sample_sort(
      tvs.begin(), tvs.end(),
      [](const voxel *a, const voxel *b) { return a->R > b->R; }, num_workers);
  std::cout << " sorting " << int(tvs.size()) << " maximal balls" << std::endl;

  ballSpace.reserve(nBalls);
  for (voxel *v : tvs) {
    ballSpace.emplace_back(v, 0);
    v->ball = &ballSpace.back();
  }

  pool.detach_blocks(0, nBalls, [&](const size_t start, const size_t end) {
    for (size_t i = start; i < end; ++i) {
      moveUphill(&ballSpace[i]);
    }
  });

  pool.wait();

  for (size_t i = 0; i < nBalls; ++i) {
    moveUphillp1(&ballSpace[i]);
  }

  pool.detach_blocks(0, nBalls, [&](const size_t start, const size_t end) {
    for (size_t i = start; i < end; ++i) {
      moveUphill(&ballSpace[i]);
    }
  });

  pool.wait();

  std::cout << " creating ball hierarchy:";
  std::cout.flush();

  for (size_t i = 0; i < nBalls; ++i) {
    findBoss(&ballSpace[i]); // 直接使用索引访问
    if (i % 100000 == 0)     // 每 100,000 个球打印一次进度
    {
      std::cout << "\r   ball: " << i;
      std::cout.flush();
    }
  }
  std::cout << "\r   ball: " << nBalls << std::endl; // 最后打印总数
}