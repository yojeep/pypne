
#include "medialSurf.h"
#include "edt.hpp"
#include "globals.h"
#include "inputData.h"
#include "typses.h"
#include <atomic>
#include <boost/multi_array.hpp>
#include <boost/sort/sort.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

medialSurface::medialSurface(
    inputDataNE &cfg) //, double vmvLimRelF, double crossAreaf
    : cg_(cfg), segs_(cfg.segs_), ToBeAssigned(0) {
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

void medialSurface::buildvoxelspace() { ///  Build voxelspace -- memory for
                                        ///  void/active voxels
  std::cout << "\nProcessing " << cg_._rockTypes[0].name
            << " voxels:" << std::endl;
  std::cout << " Creating " << nVxls << " voxels with index: " << int(0);
  std::cout.flush();

  vxlSpace.resize(nVxls);

  std::vector<voxel>::iterator p = vxlSpace.begin();
  const std::vector<voxel>::iterator vxlBegin = vxlSpace.begin();
  iZ.resize(nz);
  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      const segments &s = segs_[iz * ny + iy];
      for (int ix = 0; ix < s.cnt; ++ix) {
        if (s.s[ix].value == 0)
          for (int i = s.s[ix].start; i < s.s[ix + 1].start; ++i) {
            p->i = i;
            p->j = iy;
            p->k = iz;
            size_t currentIndex = p - vxlBegin;
            iZ[iz].push_back(currentIndex);
            ++p;
          }
      }
    }
  }

  if (nVxls != size_t(p - vxlSpace.begin()))
    std::cout << "\n Error created " << size_t(p - vxlSpace.begin())
              << " voxels " << std::endl;

  std::cout << std::endl;

  /// Link voxels to segments
  p = vxlSpace.begin();
  int total_iterations = nz * ny;
  for (int izy = 0; izy < total_iterations; ++izy) {
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
{ /// Remove maximal-balls, leave one in each adjacent voxesl. This saves time
  /// when sorting in paradoxremoveincludedballI()
  if (!nVxls) {
    return;
  }
  std::cout << " pre-remove included balls: out of " << vxlSpace.size();
  std::cout.flush();

  auto &pool = GlobalThreadPool::get();
  const int kk_steps = (nz + 1) / 2;
  const int jj_steps = (ny + 1) / 2;
  const int total_iterations = kk_steps * jj_steps;
  std::atomic<uint64_t> ndel(0);
  pool.detach_blocks(
      0, total_iterations, [&](const size_t start, const size_t end) {
        int local_ndel = 0;
        for (size_t iter = start; iter < end; ++iter) {
          int kk_idx = iter / jj_steps;
          int jj_idx = iter % jj_steps;
          int kk = kk_idx * 2;
          int jj = jj_idx * 2;
          const segments &s = segs_[kk * ny + jj];
          for (int ix = 0; ix < s.cnt; ++ix) {
            if (s.s[ix].value == 0) {
              int ii_start = s.s[ix].start;
              int ii_end = s.s[ix + 1].start;
              for (int ii = ii_start; ii < ii_end; ii += 2) {
                voxel *adjacent[8] = {
                    vxl(ii, jj, kk),         vxl(ii + 1, jj, kk),
                    vxl(ii, jj + 1, kk),     vxl(ii + 1, jj + 1, kk),
                    vxl(ii, jj, kk + 1),     vxl(ii + 1, jj, kk + 1),
                    vxl(ii, jj + 1, kk + 1), vxl(ii + 1, jj + 1, kk + 1)};

                // 查找 R 最大的体素
                voxel *max_voxel = nullptr;
                float max_R = 0.0f;
                for (int c = 0; c < 8; ++c) {
                  voxel *adj = adjacent[c];
                  if (adj != nullptr && adj->ball != nullptr &&
                      adj->R > max_R) {
                    max_voxel = adj;
                    max_R = adj->R;
                  }
                }

                // 清除非最大体素的 ball 指针
                if (max_voxel != nullptr) {
                  for (int c = 0; c < 8; ++c) {
                    voxel *adj = adjacent[c];
                    if (adj != nullptr && adj->ball != nullptr &&
                        adj != max_voxel) {
                      adj->ball = nullptr;
                      ++local_ndel;
                    }
                  }
                }
              }
            }
          }
        }
        ndel.fetch_add(local_ndel, std::memory_order_relaxed);
      });

  pool.wait();
  ndel = ndel.load();
  nBalls -= ndel;
  std::cout << ",   removed = " << ndel << " remained = " << nBalls
            << std::endl;
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
  // sort(tvs.begin(), tvs.end(), metaballcomparer());
  // sort(tvs.begin(), tvs.end(), [](const voxel *a, const voxel *b)
  // 	 { return a->R > b->R; });
  boost::sort::sample_sort(
      tvs.begin(), tvs.end(),
      [](const voxel *a, const voxel *b) { return a->R > b->R; }, num_workers);

  std::cout << " remove included balls:";
  std::cout.flush();

  size_t ndel = 0;
  std::vector<voxel *>::iterator vpp = tvs.begin(), end = tvs.end();
  while (vpp < end) {

    voxel *vi = *vpp;
    if (!vi->ball) {
      ++vpp;
      continue;
    }

    const int x = vi->i;
    const int y = vi->j;
    const int z = vi->k;
    const float ri = vi->R;
    const float ripinc = ri + 0.55; //.+RPreDelete
    const float mbmbDist = _RCorsnf * ri + _RCorsn;
    const float ripincsqr = ripinc * ripinc;
    int ex, ey, ez;
    ex = ripinc;
    for (int a = -ex; a <= ex; ++a) {
      const float asqr = a * a;
      float arg_ey = ripincsqr - asqr;
      if (arg_ey < 0)
        continue;
      ey = std::sqrtf(arg_ey);
      for (int b = -ey; b <= ey; ++b) {
        const float bsqr = b * b;
        float arg_ez = ripincsqr - asqr - bsqr;
        if (arg_ez < 0)
          continue;
        ez = std::sqrtf(arg_ez); // sqrts(r2i)+1-a-b;
        for (int c = -ez; c <= ez; ++c) {
          voxel *vj = vxl(x + a, y + b, z + c);
          if ((vj != nullptr) && (vj->ball) && (vj != vi)) {
            const float rj = vj->R;
            if (rj <= ri) {
              const float D = std::sqrtf(asqr + bsqr + c * c);
              if (D < mbmbDist || (D + rj < ripinc + _MSNoise)) {
                vj->ball = nullptr;
                ++ndel;
              }
            }
          }
        }
      }
    }

    ++vpp;
    if ((vpp - tvs.begin()) % 10000 == 0)
      std::cout << "\r  remove = " << ndel;
  }
  std::cout << "\r  removed = " << ndel << " remained = " << tvs.size() - ndel
            << " balls" << std::endl;
  nBalls -= ndel;
  return;
}

void medialSurface::moveUphill(medialBall *b_i) // const
{ /// Refines the maximal-sphere location and radius,

  const voxel *vi = vxl(b_i->fi, b_i->fj, b_i->fk);
  dbl3 disp(0., 0., 0.);
  {
    const voxel *vjm = vxl(vi->i - 1, vi->j, vi->k);
    const voxel *vjp = vxl(vi->i + 1, vi->j, vi->k);
    if (vjm && vjp) {
      float gp = vjp->R - vi->R;
      float gm = vi->R - vjm->R;
      if (std::abs(gp - gm) > 0.01)
        disp.x = std::max(-0.49, std::min(0.49, -0.5 * (gp + gm) / (gp - gm)));
    }
  }
  {
    const voxel *vjm = vxl(vi->i, vi->j - 1, vi->k);
    const voxel *vjp = vxl(vi->i, vi->j + 1, vi->k);
    if (vjm && vjp) {
      float gp = vjp->R - vi->R;
      float gm = vi->R - vjm->R;
      if (std::abs(gp - gm) > 0.01)
        disp.y = std::max(-0.49, std::min(0.49, -0.5 * (gp + gm) / (gp - gm)));
    }
  }
  {
    const voxel *vjm = vxl(vi->i, vi->j, vi->k - 1);
    const voxel *vjp = vxl(vi->i, vi->j, vi->k + 1);
    if (vjm && vjp) {
      float gp = vjp->R - vi->R;
      float gm = vi->R - vjm->R;
      if (std::abs(gp - gm) > 0.01)
        disp.z = std::max(-0.49, std::min(0.49, -0.5 * (gp + gm) / (gp - gm)));
    }
  }
  if (b_i != b_i->boss) {
    dbl3 BosKidVec = *b_i - *(b_i->boss);
    disp -=
        0.95 * ((BosKidVec & disp) / (magSqr(BosKidVec) + 1e-12)) * BosKidVec;
  }
  b_i->fi = vi->i - _mp5 + disp.x;
  b_i->fj = vi->j - _mp5 + disp.y;
  b_i->fk = vi->k - _mp5 + disp.z;
  b_i->R = vi->R + 0.95 * mag(disp);
}

void medialSurface::moveUphillp1(medialBall *bi) // const
{ /// Refines the maximal-sphere location, by moving it uphil the gradient of
  /// the distance-map, potentially relocates to new voxels

  const voxel *vi = vxl(bi->fi, bi->fj, bi->fk);
  dbl3 disp(0., 0., 0.), grad(0., 0., 0.);

  {
    const voxel *vjm = vxl(vi->i - 1, vi->j, vi->k);
    const voxel *vjp = vxl(vi->i + 1, vi->j, vi->k);
    if (vjm && vjp) {
      float gp = vjp->R - vi->R;
      float gm = vi->R - vjm->R;
      grad.x = 0.5 * (gp + gm);
      if (std::abs(gp - gm) > 0.01)
        disp.x = std::max(-0.59, std::min(0.59, -0.5 * (gp + gm) / (gp - gm)));
    }
  }
  {
    const voxel *vjm = vxl(vi->i, vi->j - 1, vi->k);
    const voxel *vjp = vxl(vi->i, vi->j + 1, vi->k);
    if (vjm && vjp) {
      float gp = vjp->R - vi->R;
      float gm = vi->R - vjm->R;
      grad.y = 0.5 * (gp + gm);
      if (std::abs(gp - gm) > 0.01)
        disp.y = std::max(-0.59, std::min(0.59, -0.5 * (gp + gm) / (gp - gm)));
    }
  }
  {
    const voxel *vjm = vxl(vi->i, vi->j, vi->k - 1);
    const voxel *vjp = vxl(vi->i, vi->j, vi->k + 1);
    if (vjm && vjp) {
      float gp = vjp->R - vi->R;
      float gm = vi->R - vjm->R;
      grad.z = 0.5 * (gp + gm);
      if (std::abs(gp - gm) > 0.01)
        disp.z = std::max(-0.59, std::min(0.59, -0.5 * (gp + gm) / (gp - gm)));
    }
  }
  disp += 1.4 * grad;

  if (bi != bi->boss) {
    dbl3 BosKidVec = *bi - *(bi->boss);
    disp -=
        0.5 * ((BosKidVec & disp) / (magSqr(BosKidVec) + 1e-12)) * BosKidVec;
  }
  disp /= (0.55 * mag(disp) + 0.05);

  voxel *vxlj = vxl(bi->fi + disp[0], bi->fj + disp[1], bi->fk + disp[2]);
  if (vxlj && vi != vxlj && vxlj->R > vi->R && vxlj->ball == nullptr) {
    bi->fi = vxlj->i - _mp5;
    bi->fj = vxlj->j - _mp5;
    bi->fk = vxlj->k - _mp5;
    bi->R = vxlj->R;
    bi->vxl->ball = nullptr;
    bi->vxl = vxlj;
    vxlj->ball = bi;
    //++nrelocations;
  }

  // cout<<nrelocations<<" relocations "<<endl;
}

void makeFriend(medialBall *vi, medialBall *vj) {
  if (vj->R > vi->R) {
    medialBall *tmp = vi;
    vi = vj;
    vj = tmp;
  }
  if ((vi->R < 1.5 * vj->R) && (!vi->isNei(vj)) && (!vi->inParents(vj)) &&
      !vj->inParents(vi)) {
    vi->addNei(vj);
    vj->addNei(vi);
  }
}

/*inline double cosAngleWithBossPerD(const medialBall* a, const medialBall* b) {
        dbl3 v1=*a-*b;
        if (b==b->boss) return 1./(mag(v1));
        dbl3 v2=*b-*(b->boss);
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
  const voxel *middlevxl = vxl(wsinv * (vi->fi * rjSqr + vj->fi * riSqr),
                               wsinv * (vi->fj * rjSqr + vj->fj * riSqr),
                               wsinv * (vi->fk * rjSqr + vj->fk * riSqr));

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

  const float x = vi->fi, y = vi->fj, z = vi->fk;
  const float ripp = vi->R * 0.6 + 2. * _MSNoise + 2.;
  const float ex = x + ripp;
  const float ripp2 = ripp * ripp;
  for (float xpa = 2.f * x - ex; xpa <= ex; xpa += 1.0f) {
    float xpa2 = (xpa - x) * (xpa - x);
    const float remain_y = ripp2 - xpa2;
    if (remain_y <= 0.f)
      continue;
    float ey = y + sqrtf(remain_y);
    for (float ypb = 2.f * y - ey; ypb <= ey; ypb += 1.0f) {
      float ypb2 = (ypb - y) * (ypb - y);
      float remain_z = remain_y - ypb2;
      if (remain_z <= 0.f)
        continue;
      float ez = z + sqrtf(remain_z);
      for (float zpc = 2.f * z - ez; zpc <= ez; zpc += 1.0f) {
        voxel *vj = this->vxl(xpa, ypb, zpc);
        if ((vj != nullptr) && vj->ball &&
            (vi !=
             vj->ball)) { //--------------------------------------------------------
          competeForParent(vi, vj->ball);
        } //--------------------------------------------------------
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
    return;
  }

  voxelImage vxls = segToVxlMesh(*this);

  auto &pool = GlobalThreadPool::get();
  uint64_t nvoxels = vxls.size();

  auto src = vxls.data_.data();
  auto inv_labels = std::make_unique_for_overwrite<bool[]>(nvoxels);
  pool.detach_blocks(0, nvoxels,
                     [&](const size_t start_idx, const size_t end_idx) {
                       for (size_t idx = start_idx; idx < end_idx; ++idx)
                         inv_labels[idx] = (src[idx] == 0);
                     });
  pool.wait();

  auto dt2 = std::make_unique_for_overwrite<float[]>(nvoxels);

  edt::binary_edtsq<bool>(inv_labels.get(), static_cast<int64_t>(vxls.nx()),
                          static_cast<int64_t>(vxls.ny()),
                          static_cast<int64_t>(vxls.nz()), 1.0f, 1.0f, 1.0f,
                          false, num_workers, dt2.get());

  const double clipROutyz = _clipROutyz;
  const double clipROutx = _clipROutx;

  std::atomic<double> totalR(0.0);

  pool.detach_blocks(
      0, nVxls, [&](const size_t start_idx, const size_t end_idx) {
        double local_R = 0.0;
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
          voxel &vit = vxlSpace[idx];
          const int x = vit.i, y = vit.j, z = vit.k;

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

          vit.R = limit;
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

  std::vector<float> delRrr(vxlSpace.size(), 0.0f);
  (std::cout << "*").flush();

  auto &pool = GlobalThreadPool::get();
  const size_t total_tasks = nz * ny;
  pool.detach_blocks(0, total_tasks, [&](const size_t start, const size_t end) {
    for (size_t iter = start; iter < end; ++iter) {
      int k = iter / ny;
      int j = iter % ny;
      const segments &s = cg_.segs_[iter];
      for (int ix = 0; ix < s.cnt; ++ix)
        if (s.s[ix].value == 0) {
          segment &seg = s.s[ix];
          for (int i = seg.start; i < s.s[ix + 1].start; ++i) {
            double sumR = 0.;
            int counter = 0;
            for (int kk = std::max(k - 1, 0); kk < std::min(k + 2, nz); ++kk)
              for (int jj = std::max(j - 1, 0); jj < std::min(j + 2, ny);
                   ++jj) {
                int ii = std::max(i - 1, 0);
                const segment *segbc = cg_.segptr(ii, jj, kk);
                if (segbc->value != 0 && (segbc + 1)->value == 0) {
                  ++segbc;
                  ii = segbc->start;
                }
                if (segbc->value == 0) {
                  int ii2 = std::min((segbc + 1)->start, i + 2);
                  voxel *vxlj = segbc->segV + (ii - segbc->start);
                  for (; ii < ii2; ++ii) {
                    sumR += vxlj->R;
                    ++vxlj;
                    counter += 1;
                  }
                }
              }

            delRrr[seg.segV + (i - seg.start) - (&vxlSpace[0])] =
                4. * sumR / (3 * counter + 27) - seg.segV[i - seg.start].R;
          }
        }
    }
  });
  pool.wait();
  (std::cout << "*").flush();

  pool.detach_blocks(0, total_tasks, [&](const size_t start, const size_t end) {
    for (size_t iter = start; iter < end; ++iter) {
      int k = iter / ny;
      int j = iter % ny;
      const segments &s = cg_.segs_[iter];
      for (int ix = 0; ix < s.cnt; ++ix)
        if (s.s[ix].value == 0) {
          segment &seg = s.s[ix];
          for (int i = seg.start; i < s.s[ix + 1].start; ++i) {
            double sumDelR = 0.;
            int counter = 0;
            for (int kk = std::max(k - 1, 0); kk < std::min(k + 2, nz); ++kk)
              for (int jj = std::max(j - 1, 0); jj < std::min(j + 2, ny);
                   ++jj) {
                int ii = std::max(i - 1, 0);
                const segment *segbc = cg_.segptr(ii, jj, kk);
                if (segbc->value != 0 && (segbc + 1)->value == 0) {
                  ++segbc;
                  ii = segbc->start;
                }
                if (segbc->value == 0) {
                  int ii2 = std::min((segbc + 1)->start, i + 2);
                  voxel *vxlj = segbc->segV + (ii - segbc->start);
                  for (; ii < ii2; ++ii) {
                    sumDelR += delRrr[vxlj - (&vxlSpace[0])];
                    ++vxlj;
                    ++counter;
                  }
                }
              }

            seg.segV[i - seg.start].R += std::min(
                std::max(
                    0.02 *
                        (delRrr[seg.segV + (i - seg.start) - (&vxlSpace[0])] -
                         0.99 * 2. * sumDelR / (1 * counter + 27)),
                    -0.005),
                0.01);
          }
        }
    }
  });
  pool.wait();
  (std::cout << "*").flush();

  /// Finally, report max distance map, to confirm that distance map is not
  /// changed too much
  float maxrrr = 0.0f;
  const size_t total = vxlSpace.size();
  if (total == 0) {
    std::cout << " maxrrr 0." << std::endl;
    return;
  }
  auto maxrrr_future = pool.submit_blocks(
      0, total, [&](const size_t start, const size_t end) -> float {
        float local_max = 0.0f;
        for (size_t i = start; i < end; ++i) {
          local_max = std::max(local_max, vxlSpace[i].R);
        }
        return local_max;
      });
  // 等待所有任务完成，并合并结果
  for (auto &future : maxrrr_future) {
    maxrrr = std::max(maxrrr, future.get());
  }
  std::cout << " maxrrr " << maxrrr << std::endl;
}

void medialSurface::
    createBallsAndHierarchy() { /// Create distance map, maximal-spheres, and
                                /// their hirarchy (medial-surface connectivity)

  //	mediaAxes medAxis(*this, minRP, clipOutSideBallFraction,
  // clipOutSideBallFraction*0.5+0.);

  buildvoxelspace();

  calc_distmaps();

  for (int i = 0; i < _nRSmoothing; ++i)
    smoothRadius();

  std::atomic<size_t> _nBalls(0);
  std::atomic<double> rBalls(0.0);

  auto &pool = GlobalThreadPool::get();
  const size_t total_iterations = vxlSpace.size();
  pool.detach_blocks(0, total_iterations,
                     [&](const size_t start, const size_t end) {
                       size_t local_nBalls = 0;
                       double local_rBalls = 0.0;

                       for (size_t i = start; i < end; ++i) {
                         voxel &v = vxlSpace[i];
                         if (v.R >= _minRp) {
                           v.ball = &ToBeAssigned;
                           ++local_nBalls;
                           local_rBalls += v.R;
                         } else {
                           v.ball = nullptr;
                         }
                       }
                       _nBalls.fetch_add(local_nBalls);
                       rBalls.fetch_add(local_rBalls);
                     });
  pool.wait();
  nBalls = _nBalls.load();

  std::cout << "\n  number of potential maximal spheres: " << nBalls
            << ",  average radius = " << rBalls.load() / nBalls << std::endl;

  paradox_pre_removeincludedballI();

  paradoxremoveincludedballI();

  std::cout << " collecting maximal balls out of " << nBalls << std::endl;

  std::vector<voxel *> tvs;
  tvs.reserve(nBalls);
  {
    for (voxel &v : vxlSpace) {
      if (v.ball) {
        tvs.push_back(&v);
      }
    }
    // sort(tvs.begin(), tvs.end(), [](const voxel *a, const voxel *b)
    // 	 { return a->R > b->R; });
    boost::sort::sample_sort(
        tvs.begin(), tvs.end(),
        [](const voxel *a, const voxel *b) { return a->R > b->R; },
        num_workers);
    std::cout << " sorting " << int(tvs.size()) << " maximal balls"
              << std::endl;
  }

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