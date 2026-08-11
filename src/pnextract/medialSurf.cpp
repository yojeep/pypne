#include <algorithm>
#include <boost/sort/sort.hpp>
#include <iostream>
#include <vector>
#include "medialSurf.hpp"
#include "edt_float_bv.hpp"
#include "globals.hpp"
#include "indexunraveler.hpp"
#include "types.hpp"

void medialSurface::setDefaults(float avgR) {
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
  if (cfg.minRp.has_value())
    _minRp = cfg.minRp.value();
  else if (avgR > 0) {
    std::cout << "keyword \"minRp\" not found, default value "
              << std::abs(_minRp) << " will be used\n"
              << std::endl;
  }
  // before distance map
  _clipROutx = 0.05;
  _clipROutyz = 0.98;
  // after distance map
  _midRf = 0.7;
  _MSNoise = 1. * std::abs(_minRp) + 1.;
  _lenNf = 0.6;
  _vmvRadRelNf = 1.1;
  _nRSmoothing = 3;
  _RCorsnf = 0.15;
  _RCorsn = std::abs(_minRp);

  if (cfg.nBP6 == 6)
    _clipROutyz = _clipROutx;

  if (cfg.clipROutx.has_value())
    _clipROutx = cfg.clipROutx.value();
  if (cfg.clipROutyz.has_value())
    _clipROutyz = cfg.clipROutyz.value();
  if (cfg.midRf.has_value())
    _midRf = cfg.midRf.value();
  if (cfg.MSNoise.has_value())
    _MSNoise = cfg.MSNoise.value();
  if (cfg.lenNf.has_value())
    _lenNf = cfg.lenNf.value();
  if (cfg.vmvRadRelNf.has_value())
    _vmvRadRelNf = cfg.vmvRadRelNf.value();
  if (cfg.nRSmoothing.has_value())
    _nRSmoothing = cfg.nRSmoothing.value();
  if (cfg.RCorsnf.has_value())
    _RCorsnf = cfg.RCorsnf.value();
  if (cfg.RCorsn.has_value())
    _RCorsn = cfg.RCorsn.value();

  if (_minRp < 0.)
    return;

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

void medialSurface::buildvoxelspace() { ///  Build voxelspace -- memory for
                                        ///  void/active voxels

  std::cout << " Creating " << nVxls << " voxels with index: " << int(1);
  std::cout.flush();
  for (int z = 0; z < nz; ++z)
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x)
        if (bool_img(z, y, x)) {
          vxlMap.vec_data_.emplace_back(z, y, x, 0.0f);
        }
}

void medialSurface::paradox_pre_removeincludedballI() // to remove the included
                                                      // maximal bals
{ /// Remove maximal-balls, leave one in each adjacent voxel. This saves time
  /// when sorting in paradoxremoveincludedballI()

  auto &pool = GlobalThreadPool::get();
  const size_t z_steps = (nz + 1) / 2;
  const size_t y_steps = (ny + 1) / 2;
  const size_t x_steps = (nx + 1) / 2;
  const size_t total_iterations = z_steps * y_steps * x_steps;

  IndexUnraveler unraveler({z_steps, y_steps, x_steps});

  std::atomic<uint64_t> _nBalls(0);
  pool.detach_blocks(
      0, total_iterations, [&](const size_t start, const size_t end) noexcept {
        std::array<int, 3> coord{};
        for (size_t idx = start; idx < end; ++idx) {
          unraveler.unravel(idx, coord);
          int zo = coord[0] * 2, yo = coord[1] * 2, xo = coord[2] * 2;
          float max_r = -1.0f;
          voxel *max_voxel = nullptr;
          for (int zi = 0; zi < 2; ++zi)
            for (int yi = 0; yi < 2; ++yi)
              for (int xi = 0; xi < 2; ++xi) {
                int z = zo + zi, y = yo + yi, x = xo + xi;
                if (!bool_img.isinside(z, y, x) || !bool_img(z, y, x))
                  continue;
                voxel &v = vxlMap(z, y, x);
                float r = v.R;
                if (r <= max_r && r <= _minRp)
                  continue;
                max_r = r;
                max_voxel = &v;
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
  std::cout << " sorting... ";
  std::vector<voxel *> tvs;
  tvs.reserve(nBalls);

  for (voxel &v : vxlMap.vec_data_) {
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

    float ri = vi->R;
    float r_sphere = ri + 0.55f;
    float mbmbDist = _RCorsnf * ri + _RCorsn;
    int iz = vi->z, iy = vi->y, ix = vi->x;

    for_each_voxel_in_sphere(
        iz, iy, ix, r_sphere, nz, ny, nx, [&](int z, int y, int x) noexcept {
          if (!bool_img(z, y, x))
            return;
          voxel &vj = vxlMap(z, y, x);
          if (!vj.ball || vi == &vj)
            return;
          float rj = vj.R;
          if (rj <= ri) {
            float dz = z - iz, dy = y - iy, dx = x - ix;
            float D = norm(dz, dy, dx);
            if (D >= mbmbDist && (D + rj >= r_sphere + _MSNoise))
              return;
            vj.ball = nullptr;
            ++ndel;
          }
        });

    if ((idx % 10000) == 0)
      std::cout << "\r  remove = " << ndel;
  }
  std::cout << "\r  removed = " << ndel << " remained = " << tvs.size() - ndel
            << " balls" << std::endl;
  nBalls -= ndel;
  return;
}

void medialSurface::moveUphill(medialBall &bi) {
  const voxel &vi = vxlMap(bi.iz(), bi.iy(), bi.ix());
  const int zm = vi.z, ym = vi.y, xm = vi.x;

  const std::array<int, 3> center = {xm, ym, zm};
  const std::array<int, 3> dims = {nx, ny, nz};

  Vector3f64 disp{0.0f, 0.0f, 0.0f};

  for (int a = 0; a < 3; ++a) {
    const int cm = center[a];
    const int dimSize = dims[a];

    if (cm < 1 || cm > dimSize - 2)
      continue;

    const auto &off = Constants::NEIGHBOR_OFFSETS_3[a];
    const voxel &vjn = vxlMap(zm - off[0], ym - off[1], xm - off[2]); // neg
    const voxel &vjp = vxlMap(zm + off[0], ym + off[1], xm + off[2]); // pos

    if (!vjn || !vjp)
      continue;

    const float gp = vjp.R - vi.R;
    const float gm = vi.R - vjn.R;
    const float denom = gp - gm;

    if (std::abs(denom) < 0.01f)
      continue;

    disp[a] = std::clamp(-0.5f * (gp + gm) / denom, -0.49f, 0.49f);
  }

  if (bi.boss != &bi) {
    const Vector3f64 BosKidVec = {bi.fx - bi.boss->fx, bi.fy - bi.boss->fy,
                                  bi.fz - bi.boss->fz};
    const float sqNorm = BosKidVec.squaredNorm() + Constants::F32_EPS;
    disp -= 0.95f * (BosKidVec.dot(disp) / sqNorm) * BosKidVec;
  }

  bi.fx = vi.x + _0p5 + disp[0];
  bi.fy = vi.y + _0p5 + disp[1];
  bi.fz = vi.z + _0p5 + disp[2];
  bi.R = vi.R + 0.95f * disp.norm();
}

void medialSurface::moveUphillp1(medialBall &bi) {
  /// Refines the maximal-sphere location by moving it uphill on the
  /// distance-map gradient, potentially relocating to a new voxel.
  const voxel &vi = vxlMap(bi.iz(), bi.iy(), bi.ix());
  const int zm = vi.z, ym = vi.y, xm = vi.x;

  const std::array<int, 3> center = {xm, ym, zm};
  const std::array<int, 3> dims = {nx, ny, nz};

  Vector3f64 disp{0.0f, 0.0f, 0.0f};
  Vector3f64 grad{0.0f, 0.0f, 0.0f};

  for (int a = 0; a < 3; ++a) {
    const int cm = center[a];
    const int dimSize = dims[a];

    if (cm < 1 || cm > dimSize - 2)
      continue;

    const auto &off = Constants::NEIGHBOR_OFFSETS_3[a];
    const voxel &vjm = vxlMap(zm - off[0], ym - off[1], xm - off[2]);
    const voxel &vjp = vxlMap(zm + off[0], ym + off[1], xm + off[2]);

    const float gp = vjp.R - vi.R;
    const float gm = vi.R - vjm.R;
    const float sum = gp + gm;
    const float denom = gp - gm;

    grad[a] = 0.5f * sum;

    if (std::abs(denom) < 0.01f)
      continue;
    disp[a] = std::clamp(-0.5f * sum / denom, -0.59f, 0.59f);
  }

  disp += 1.4f * grad;

  if (bi.boss != &bi) {
    const Vector3f64 BosKidVec = {bi.fx - bi.boss->fx, bi.fy - bi.boss->fy,
                                  bi.fz - bi.boss->fz};
    const float sqNorm = BosKidVec.squaredNorm() + Constants::F32_EPS;
    disp -= 0.5f * (BosKidVec.dot(disp) / sqNorm) * BosKidVec;
  }

  disp /= 0.55f * disp.norm() + 0.05f;

  const int ix_vj = static_cast<int>(bi.fx + disp[0]);
  const int iy_vj = static_cast<int>(bi.fy + disp[1]);
  const int iz_vj = static_cast<int>(bi.fz + disp[2]);

  if (!vxlMap.isinside(iz_vj, iy_vj, ix_vj)) {
    return;
  }
  voxel &vj = vxlMap(iz_vj, iy_vj, ix_vj);
  if (&vi == &vj || vj.R <= vi.R || vj.ball != nullptr)
    return;
  bi.vxl->ball = nullptr;
  bi.vxl = &vj;
  bi.fx = vj.x + _0p5;
  bi.fy = vj.y + _0p5;
  bi.fz = vj.z + _0p5;
  bi.R = vj.R;
  vj.ball = &bi;
}

void medialSurface::competeForParent(medialBall *bi, medialBall *bj) {

  auto ratioBoss = [](float bossR, float selfR, float d, float n) noexcept {
    return (bossR - selfR + 2.0f * n) / (d + 0.25f);
  };
  auto ratioPeer = [](float otherR, float selfR, float d, float n) noexcept {
    return (otherR - selfR + 2.0f * n + 0.01f) / (d + 0.2f);
  };
  auto ratioNear = [](float bossR, float selfR, float d, float n,
                      float off) noexcept {
    return (bossR - selfR + 2.0f * n) / (d + off);
  };

  const float noise = _MSNoise;
  const float ri = bi->R;
  const float rj = bj->R;
  const float riSq = ri * ri;
  const float rjSq = rj * rj;
  const float dSq = distSq(bi, bj);
  const float dSqrt = std::sqrt(dSq);

  // Calculate weighted middle voxel
  const float wsinv = 1.0f / (riSq + rjSq);
  int iz_middle = wsinv * (bi->fz * rjSq + bj->fz * riSq),
      iy_middle = wsinv * (bi->fy * rjSq + bj->fy * riSq),
      ix_middle = wsinv * (bi->fx * rjSq + bj->fx * riSq);

  if (!vxlMap.isinside(iz_middle, iy_middle, ix_middle))
    return;
  const voxel &middlevxl = vxlMap(iz_middle, iy_middle, ix_middle);
  if (!middlevxl || middlevxl.R < std::min(ri, rj) * _midRf - 0.5f ||
      1.01 * dSqrt > ri + rj + 1.0f + noise)
    return;

  // Initialize Boss assignment
  if (bj->boss == bj && bi->mastrSphere() != bj) {
    if (ri >= rj)
      bj->boss = bi;
    else if (bi->boss->R <= rj)
      bi->boss = bj;
    else if (ri >= rj - noise && ri * _vmvRadRelNf + noise >= rj)
      bj->boss = bi;
  } else if (bi->boss == bi && bj->mastrSphere() != bi) {
    if (rj >= ri)
      bi->boss = bj;
    else if (bj->boss->R <= ri)
      bj->boss = bi;
    else if (rj >= ri - noise && rj * _vmvRadRelNf + noise >= ri)
      bi->boss = bj;
  }

  medialBall *mbi = bi->mastrSphere();
  medialBall *mbj = bj->mastrSphere();

  if (mbi == bj || mbj == bi) // If already boss and slave, return
    return;

  if (mbi == mbj) {
    {
      const int leveli = bi->level();
      const int levelj = bj->level();
      const int levelDiff = leveli - levelj;

      const float distBiBoss = dist(bi->boss, bi);
      const float distBjBoss = dist(bj->boss, bj);
      const float distBiBj = dSqrt;

      if (levelDiff < -1 && ratioBoss(bj->boss->R, bj->R, distBjBoss, noise) <
                                ratioPeer(bi->R, bj->R, distBiBj, noise))
        bj->boss = bi;
      else if (levelDiff > 1 &&
               ratioBoss(bi->boss->R, bi->R, distBiBoss, noise) <
                   ratioPeer(bj->R, bi->R, distBiBj, noise))
        bi->boss = bj;
      else if (levelDiff > 0 &&
               ratioNear(bi->boss->R, bi->R, distBiBoss, noise, 1.2f) <
                   ratioNear(bj->R, bi->R, distBiBj, noise, 1.3f) &&
               !bj->inParents(bi))
        bi->boss = bj;
      else if (levelDiff < 0 &&
               ratioNear(bj->boss->R, bj->R, distBjBoss, noise, 1.2f) <
                   ratioNear(bi->R, bj->R, distBiBj, noise, 1.3f) &&
               !bi->inParents(bj))
        bj->boss = bi;
    }
    if (bi->mastrSphere() != bj->mastrSphere()) {
      std::cout << "sdsdsds" << std::endl;
      exit(-1);
    }
  } else { // mvi != mvj: Merge different master balls
    const float avgR = 0.5f * (mbi->R + mbj->R);
    if (distSq(mbi, mbj) < _lenNf * sq(avgR + 2.0f * noise)) {
      // Ensure mvi is the larger one
      if (mbi->R < mbj->R) {
        std::swap(bi, bj);
        std::swap(mbi, mbj);
      }
      // Merge small tree nodes to larger tree
      if (mbj->R < _vmvRadRelNf * bj->R + noise &&
          mbj->R < _vmvRadRelNf * bi->R + noise &&
          mbj->R < _vmvRadRelNf * bi->boss->R + noise) {
        while (bj->boss != bj && mbj->R < _vmvRadRelNf * bj->boss->R + noise) {
          medialBall *tmp = bj->boss;
          bj->boss = bi;
          bi = bj;
          bj = tmp;
        }
        if (bj->boss == bj && bi->mastrSphere() != bj) {
          bj->boss = bi;
        }
      }
    }

    // Final level balance
    if (bi != bj->boss) {
      mbi = bi->mastrSphere();
      mbj = bj->mastrSphere();

      int leveli = bi->level();
      int levelj = bj->level();
      const float distAvg = dist(mbi, mbj) + 0.5f * noise;

      // Symmetric level-balancing step extracted as lambda
      auto balanceStep = [&](medialBall *&ba, medialBall *&bb, medialBall *&ma,
                             medialBall *&mb, int &levA, int &levB) noexcept {
        if (levA < levB)
          return false;
        const float dMaPa = dist(ma, ba);
        const float dMbPa = dist(mb, ba);
        if ((ba->boss->R - ba->R + 0.55f * noise) / (dMaPa + distAvg) >=
            (bb->R - ba->R + 0.5f * noise) / (dMbPa + distAvg))
          return false;
        medialBall *tmp = ba->boss;
        ba->boss = bb;
        bb = ba;
        ba = tmp;
        ++levB;
        --levA;
        return true;
      };

      while (balanceStep(bi, bj, mbi, mbj, leveli, levelj)) {
      }
      while (balanceStep(bj, bi, mbj, mbi, levelj, leveli)) {
      }
    }
  }
}

void medialSurface::findBoss(medialBall *bi) {

  float zo = bi->fz, yo = bi->fy, xo = bi->fx;
  float r_sphere = bi->R * 0.6f + 2.0f * _MSNoise + 2.0f;

  for_each_voxel_in_sphere(zo, yo, xo, r_sphere, nz, ny, nx,
                           [&](int z, int y, int x) noexcept {
                             if (!bool_img(z, y, x))
                               return;
                             voxel &vj = vxlMap(z, y, x);
                             if (!vj.ball || bi == vj.ball)
                               return;
                             competeForParent(bi, vj.ball);
                           });
}

float medialSurface::calc_distmaps() // search  MBs at each voxel
{

  std::cout << " computing distance map for index " << int(0);
  std::cout.flush();

  std::unique_ptr<float[]> dt2 = edt::binary_edtsq(
      bool_img.bv_, static_cast<int64_t>(nx), static_cast<int64_t>(ny),
      static_cast<int64_t>(nz), 1.0f, 1.0f, 1.0f, false);

  std::cout << "  distance map computed" << std::endl;

  const float clipROutyz = _clipROutyz;
  const float clipROutx = _clipROutx;

  std::atomic<double> totalR(0.0);
  auto &pool = GlobalThreadPool::get();
  pool.detach_blocks(
      0, nVxls, [&](const size_t start_idx, const size_t end_idx) noexcept {
        double local_R = 0.0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          voxel &vi = vxlMap.vec_data_[idx];
          const int z = vi.z, y = vi.y, x = vi.x;

          float limit = std::sqrtf(dt2[(z * ny + y) * nx + x]) - 0.5f;

          float iSqr = std::min(y + 2.0f, ny - y + 1.0f);
          if (iSqr < limit)
            limit = std::max((1.0f - clipROutyz) * limit + clipROutyz * iSqr,
                             0.01f);

          iSqr = std::min(z + 2.0f, nz - z + 1.0f);
          if (iSqr < limit)
            limit = std::max((1.0f - clipROutyz) * limit + clipROutyz * iSqr,
                             0.01f);

          iSqr = std::min(x + 2.0f, nx - x + 1.0f);
          if (iSqr < limit)
            limit =
                std::max((1.0f - clipROutx) * limit + clipROutx * iSqr, 0.1f);

          vi.R = limit;
          local_R += limit;
        }
        totalR.fetch_add(local_R, std::memory_order_relaxed);
      });

  pool.wait();
  float rBalls = totalR.load();
  float avgR = rBalls / nVxls;

  std::cout << "\n";
  std::cout << "  average distance map = " << avgR << std::endl;

  return avgR;
}

void medialSurface::smoothRadius() {

  (std::cout << " smoothing R  ").flush();
  (std::cout << "*").flush();

  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = nVxls;

  Eigen::Tensor<float, 3, Eigen::RowMajor> delRrr(nz, ny, nx);

  pool.detach_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) noexcept {
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          voxel &vi = vxlMap.vec_data_[idx];
          const int zo = vi.z, yo = vi.y, xo = vi.x;
          float sumR{0.};
          size_t counter{0};
          for (int zi = -1; zi < 2; ++zi)
            for (int yi = -1; yi < 2; ++yi)
              for (int xi = -1; xi < 2; ++xi) {
                int z = zo + zi, y = yo + yi, x = xo + xi;
                if (!bool_img.isinside(z, y, x) || !bool_img(z, y, x))
                  continue;
                sumR += vxlMap(z, y, x).R;
                ++counter;
              }
          delRrr(zo, yo, xo) = 4. * sumR / (3 * counter + 27) - vi.R;
        }
      });
  pool.wait();
  (std::cout << "*").flush();

  pool.detach_blocks(
      0, total_iterations,
      [&](const size_t start_idx, const size_t end_idx) noexcept {
        for (size_t idx = start_idx; idx < end_idx; ++idx) {

          voxel &vi = vxlMap.vec_data_[idx];
          const int zo = vi.z, yo = vi.y, xo = vi.x;
          float sumDelR{0.};
          size_t counter{0};
          for (int zi = -1; zi < 2; ++zi)
            for (int yi = -1; yi < 2; ++yi)
              for (int xi = -1; xi < 2; ++xi) {
                int z = zo + zi, y = yo + yi, x = xo + xi;
                if (!bool_img.isinside(z, y, x) || !bool_img(z, y, x))
                  continue;
                sumDelR += delRrr(z, y, x);
                ++counter;
              }
          vi.R += std::clamp(0.02f * (delRrr(zo, yo, xo) -
                                      0.99f * 2.0f * sumDelR / (counter + 27)),
                             -0.005f, 0.01f);
        }
      });
  pool.wait();
  (std::cout << "*").flush();

  float maxrrr =
      std::max_element(vxlMap.vec_data_.begin(), vxlMap.vec_data_.end(),
                       [](const auto &a, const auto &b) { return a.R < b.R; })
          ->R;

  std::cout << " maxrrr " << maxrrr << std::endl;
}

void medialSurface::createBallsAndHierarchy() { /// Create distance map,
                                                /// maximal-spheres, and their
                                                /// hirarchy (medial-surface
                                                /// connectivity)

  buildvoxelspace();

  float avgR = calc_distmaps();

  setDefaults(avgR);

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

  for (voxel &v : vxlMap.vec_data_) {
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
  std::cout << "  number of maximal balls: " << nBalls << std::endl;

  pool.detach_blocks(0, nBalls,
                     [&](const size_t start, const size_t end) noexcept {
                       for (size_t i = start; i < end; ++i) {
                         moveUphill(ballSpace[i]);
                       }
                     });

  pool.wait();
  std::cout << "moveUphill done" << std::endl;

  for (size_t i = 0; i < nBalls; ++i) {
    moveUphillp1(ballSpace[i]);
  }
  std::cout << "moveUphillp1 done" << std::endl;

  pool.detach_blocks(0, nBalls,
                     [&](const size_t start, const size_t end) noexcept {
                       for (size_t i = start; i < end; ++i) {
                         moveUphill(ballSpace[i]);
                       }
                     });

  pool.wait();

  std::cout << " creating ball hierarchy:";
  std::cout.flush();

  for (size_t i = 0; i < nBalls; ++i) {
    medialBall *p = &ballSpace[i];
    findBoss(p);
    if (i % 100000 == 0) // 每 100,000 个球打印一次进度
    {
      std::cout << "\r   ball: " << i;
      std::cout.flush();
    }
  }
  std::cout << "\r   ball: " << nBalls << std::endl; // 最后打印总数
}