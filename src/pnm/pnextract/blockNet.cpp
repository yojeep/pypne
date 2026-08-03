#include "blockNet.h"
#include <iostream>
#include <vector>

// #include "vxlImage_manip.h"

// clock_t myTime::start = clock();

int blockNetwork::nBP6 = 2;     //   will be set in inputDataNE::init...
int blockNetwork::SideImax = 1; //=nBP6-1

void blockNetwork::createMedialSurface(
    medialSurface *&srf, inputDataNE &cfg,
    size_t startValue) { /// calls medialSurface::createBallsAndHierarchy(), ...

  {

    //  clipROutx clipROutyz  midRFrac  RMedSurfNoise  lenNf   vmvRadRelNf
    //  nRSmoothing   RCorsf   RCors

    medialSurface *medSurf = new medialSurface(cfg);
    medSurf->createBallsAndHierarchy();
    srf = (medSurf);
  }
}

void blockNetwork::CreateVElem(
    size_t
        startValue) { /// ### map pore labels from maximal spheres to the image
                      /// blockNetwork::VElems, standing for VoxelElements.
  const Eigen::TensorMap<Eigen::Tensor<uint8_t, 3, Eigen::RowMajor>>
      &binary_image = srf->binary_image;
  const std::vector<voxel> &vxlSpace = srf->vxlSpace;
  const Eigen::Tensor<voxel *, 3, Eigen::RowMajor> &vxlMap = srf->vxlMap;
  const std::vector<medialBall> &ballSpace = srf->ballSpace;

  std::cout << "\n\nCreating pore elements:" << std::endl;

  if (poreIs.empty()) {
    sides[0] = voxel(-1, cg.ny / 2, cg.nz / 2, _SideRadius); // minX
    poreNE *tmp = new poreNE();
    tmp->volumn = cg.ny * cg.nz;
    tmp->mb = new medialBall(&sides[0], -1);
    poreIs.push_back(tmp);

    sides[1] = voxel(cg.nx + 1, cg.ny / 2, cg.nz / 2, _SideRadius); // maxX
    tmp = new poreNE();
    tmp->volumn = cg.ny * cg.nz;
    tmp->mb = new medialBall(&sides[1], -1);
    poreIs.push_back(tmp);

    if (nBP6 == 6) {
      sides[2] = voxel(cg.nx / 2, -1, cg.nz / 2, _SideRadius); // minY
      tmp = new poreNE();
      tmp->volumn = cg.nx * cg.nz;
      tmp->mb = new medialBall(&sides[2], -1);
      poreIs.push_back(tmp);

      sides[3] = voxel(cg.nx / 2, cg.ny + 1, cg.nz / 2, _SideRadius); // maxY
      tmp = new poreNE();
      tmp->volumn = cg.nx * cg.nz;
      tmp->mb = new medialBall(&sides[3], -1);
      poreIs.push_back(tmp);

      sides[4] = voxel(cg.nx / 2, cg.ny / 2, -1, _SideRadius); // minZ
      tmp = new poreNE();
      tmp->volumn = cg.ny * cg.nx;
      tmp->mb = new medialBall(&sides[4], -1);
      poreIs.push_back(tmp);

      sides[5] = voxel(cg.nx / 2, cg.ny / 2, cg.nz + 1, _SideRadius); // maxZ
      tmp = new poreNE();
      tmp->volumn = cg.ny * cg.nx;
      tmp->mb = new medialBall(&sides[5], -1);
      poreIs.push_back(tmp);
    }
  }

  // VElems.reset(cg.nx + 2, cg.ny + 2, cg.nz + 2, -257);
  VElems.reset(cg.nx + 2, cg.ny + 2, cg.nz + 2, -1);
  VElems.X0Ch() = cg.VImage.X0() - cg.VImage.dx();
  VElems.dxCh() = cg.VImage.dx();
  // int nVVs = 0;
  // for (int iz = 0; iz < cg.nz; ++iz) {
  //   for (int iy = 0; iy < cg.ny; ++iy) {
  //     const segments &s = cg.segs_[iz * cg.ny + iy];
  //     for (int ix = 0; ix < s.cnt; ++ix) {
  //       int value = -1 - int(s.s[ix].value);
  //       nVVs = std::max(nVVs, value);
  //       for (int i = s.s[ix].start; i < s.s[ix + 1].start; ++i) {
  //         VElems(i + 1, iy + 1, iz + 1) = value;
  //       }
  //     }
  //   }
  // }
  // ++nVVs;
  int nVVs = 1;
  VElems.setSlice('i', 0, 0);
  VElems.setSlice('i', cg.nx + 1, 1);
  if (nBP6 == 6) {
    VElems.setSlice('j', 0, 2); // TODO: check
    VElems.setSlice('j', cg.ny + 1, 3);
    VElems.setSlice('k', 0, 4);
    VElems.setSlice('k', cg.nz + 1, 5);
  } else {
    VElems.setSlice('j', 0,
                    -1 - 2 - nVVs); // TODO: check, this means negatives can not
                                    // be used to uniquely identify pores
    VElems.setSlice('j', cg.ny + 1, -1 - 3 - nVVs);
    VElems.setSlice('k', 0, -1 - 4 - nVVs);
    VElems.setSlice('k', cg.nz + 1, -1 - 5 - nVVs);
  }

  firstPore = 2;

  {

    int uasyned = -1;

    firstPores = (poreIs.size());
    for (const auto &bi : ballSpace) {
      if (bi.boss == &bi) {
        VElems(bi.fi + 1, bi.fj + 1, bi.fk + 1) = poreIs.size();
        poreNE *tmp = new poreNE();
        tmp->mb = &bi;
        poreIs.push_back(tmp);
      }
    }
    std::cout << "\n created " << poreIs.size()
              << " pores (+boundaries)  for up to index " << 0 << std::endl;
    lastPores = (poreIs.size() - 1);

    // const int firstPoreInd = firstPores;
    // const int lastPoreInd = lastPores;
    std::cout << " mapping pores indices to image, for index " << 0 << ":  "
              << firstPores << " to " << lastPores << ",  unasigned:" << uasyned
              << std::endl;

    for (const auto &bi : ballSpace) {
      // if (!bi.boss)
      //   continue;

      medialBall *mastrSphere = bi.mastrSphere();

      int Vm =
          VElems(mastrSphere->fi + 1, mastrSphere->fj + 1, mastrSphere->fk + 1);
      // ensure(VElemV > 0 && VElemV < len(poreIs));

      const float zo = bi.fk, yo = bi.fj, xo = bi.fi;

      float r = bi.R;
      int rz2 = std::pow(std::max(r * 0.25f - 1.f, 1.001f), 2);
      float rz = std::sqrtf(rz2);
      for (float zf = std::max(zo - rz, 0.5f);
           zf < std::min(zo + rz, cg.nz - 0.5f) + 1e-6f; zf += 1.0f) {
        float ry2 = rz2 - std::pow(zf - zo, 2);
        if (ry2 <= 0.f)
          continue;
        float ry = std::sqrtf(ry2);
        for (float yf = std::max(yo - ry, 0.5f);
             yf < std::min(yo + ry, cg.ny - 0.5f) + 1e-6f; yf += 1.0f) {
          float rx2 = ry2 - std::pow(yf - yo, 2);
          if (rx2 <= 0.f)
            continue;
          float rx = std::sqrtf(rx2);
          for (float xf = std::max(xo - rx, 0.5f);
               xf < std::min(xo + rx, cg.nx - 0.5f) + 1e-6f; xf += 1.0f) {
            int z = static_cast<int>(zf);
            int y = static_cast<int>(yf);
            int x = static_cast<int>(xf);

            if (x < 0 || y < 0 || z < 0 || x >= cg.nx || y >= cg.ny ||
                z >= cg.nz || binary_image(z, y, x) != 0)
              continue;
            int zVE = z + 1;
            int yVE = y + 1;
            int xVE = x + 1;
            int Vi = VElems(xVE, yVE, zVE);
            if (Vi == uasyned)
              VElems(xVE, yVE, zVE) = Vm;
            else if (Vi != Vm) {
              voxel *vi = vxlMap(z, y, x);
              if (!vi->ball && vi->R < r) {
                const medialBall *bvi = poreIs[Vi]->mb;
                float zif = zf - bvi->fk, yif = yf - bvi->fj,
                      xif = xf - bvi->fi;
                float zmf = zf - mastrSphere->fk, ymf = yf - mastrSphere->fj,
                      xmf = xf - mastrSphere->fi;
                if (zmf * zmf + ymf * ymf + xmf * xmf <
                    zif * zif + yif * yif + xif * xif)
                  VElems(xVE, yVE, zVE) = Vm; // elem->id;
              }
            }
          }
        }
      }
    }

    voxelField<int> voxls = VElems;
    PoreGrower pg(cg, poreIs, VElems, voxls, vxlSpace, vxlMap, ballSpace,
                  firstPores, lastPores, uasyned);
    pg.grow();
  }
}

void blockNetwork::createNewThroats(medialSurface *&srf) {

  // vector<int> iThroatFaces;
  // iThroatFaces.reserve(throatIs.size()+poreIs.size()*10);
  const Eigen::TensorMap<Eigen::Tensor<uint8_t, 3, Eigen::RowMajor>>
      &binary_image = srf->binary_image;
  const std::vector<voxel> &vxlSpace = srf->vxlSpace;
  const Eigen::Tensor<voxel *, 3, Eigen::RowMajor> &vxlMap = srf->vxlMap;
  const std::vector<medialBall> &ballSpace = srf->ballSpace;
  struct NeighborOffset {
    int dz, dy, dx;
  };

  static constexpr std::array<NeighborOffset, 6> NEIGHBOR_OFFSETS = {{
      {0, 0, 1}, // i+1
      {0, 1, 0}, // j+1
      {1, 0, 0}  // k+1
  }};

  {
    std::cout << "\nlooking for connections, iFirst = " << throatIs.size()
              << " ..  ";
    std::cout.flush();

    for (voxel vi : vxlSpace) {
      int x0 = vi.i;
      int y0 = vi.j;
      int z0 = vi.k;
      int p0ID = VElems(x0, y0, z0);

      if (p0ID < 0)
        continue;
      for (const auto &offset : NEIGHBOR_OFFSETS) {
        int x1 = x0 + offset.dx;
        int y1 = y0 + offset.dy;
        int z1 = z0 + offset.dz;
        int p1ID = VElems(x1, y1, z1);
        if (p0ID == p1ID || p1ID < 0 || (p0ID + p1ID) < 2)
          continue;
        bool dir = p1ID > p0ID;
        if (!dir) {
          std::swap(p0ID, p1ID);
        }
        poreNE *p0 = poreIs[p0ID];
        poreNE *p1 = poreIs[p1ID];

        if (!p0)
          std::cout << "  ERROR p0ID" << p0ID << std::endl;
        int tIDNext = throatIs.size();
        std::pair<std::map<int, int>::iterator, bool> ret =
            p1->contacts.insert(std::pair<int, int>(p0ID, tIDNext));
        if (ret.second) {

          if (ret.second !=
              p0->contacts.insert(std::pair<int, int>(p1ID, tIDNext)).second)
            std::cout << "Errorpnm821-2 " << p0ID << " " << p1ID << std::endl;
          throatIs.push_back(new throatNE(tIDNext, p0ID, p1ID));
          // if (p1ID<2 &&p2ID<2) cout <<"Errskziw1: p1ID<2 &&p2ID<2:
          // tid: "<< tIDNext<<"  "<<i<<" "<<j<<" "<<k<<" " <<"
          // "<<p1->mb->fi<<" "<<p1->mb->fj<<" "<<p1->mb->fk<<" "<<"
          // "<<p2->mb->fi<<" "<<p2->mb->fj<<" "<<p2->mb->fk<<" "<<endl;
        }
        throatNE *trot = throatIs[ret.first->second];
        trot->CrosArea.x += (2 * dir - 1);
        // trot->C[0] += int3{{i-1,j-1,k-1}};

        if (p0ID >= firstPore) {
          ++(p0->surfaceArea);
          ++(p0->volumn);
        }
        if (p1ID >= firstPore)
          ++(p1->surfaceArea);
      }
    }

    for (int z = 1; z < VElems.nz() - 1; ++z)
      for (int y = 1; y < VElems.ny() - 1; ++y)
        for (int x = 0; x < VElems.nx() - 1; ++x) {
          int p1ID = VElems(x, y, z);
          if (p1ID >= 0) {
            int p2ID = VElems(x + 1, y, z);
            if (p1ID != p2ID) {
              if (p2ID >= 0) {
                bool dir = p2ID > p1ID;
                if (!dir) {
                  int tmp = p1ID;
                  p1ID = p2ID;
                  p2ID = tmp;
                }

                poreNE *p1 = poreIs[p1ID];
                poreNE *p2 = poreIs[p2ID];

                if (!p1)
                  std::cout << "  ERROR p1ID" << p1ID << std::endl;
                int tIDNext = throatIs.size();
                std::pair<std::map<int, int>::iterator, bool> ret =
                    p2->contacts.insert(std::pair<int, int>(p1ID, tIDNext));
                if (ret.second) {

                  if (ret.second !=
                      p1->contacts.insert(std::pair<int, int>(p2ID, tIDNext))
                          .second)
                    std::cout << "Errorpnm821-2 " << p1ID << " " << p2ID
                              << std::endl;
                  throatIs.push_back(new throatNE(tIDNext, p1ID, p2ID));
                  // if (p1ID<2 &&p2ID<2) cout <<"Errskziw1: p1ID<2 &&p2ID<2:
                  // tid: "<< tIDNext<<"  "<<i<<" "<<j<<" "<<k<<" " <<"
                  // "<<p1->mb->fi<<" "<<p1->mb->fj<<" "<<p1->mb->fk<<" "<<"
                  // "<<p2->mb->fi<<" "<<p2->mb->fj<<" "<<p2->mb->fk<<" "<<endl;
                }
                throatNE *trot = throatIs[ret.first->second];
                trot->CrosArea.x += (2 * dir - 1);
                // trot->C[0] += int3{{i-1,j-1,k-1}};
              }

              if (p1ID >= firstPore)
                ++(poreIs[p1ID]->surfaceArea);
              if (p2ID >= firstPore)
                ++(poreIs[p2ID]->surfaceArea);
            }

            if (p1ID >= firstPore) {
              ++(poreIs[p1ID]->volumn);
            }
          }
        }

    (std::cout << throatIs.size() << " .. ").flush();

    for (int z = 1; z < VElems.nz() - 1; z++) //     y >-<
      for (int y = 0; y < VElems.ny() - 1; ++y)
        for (int x = 1; x < VElems.nx() - 1; ++x) {
          int p1ID = VElems(x, y, z);
          if (p1ID >= 0) {
            int p2ID = VElems(x, y + 1, z);
            if (p1ID != p2ID) {
              if (p2ID >= 0) {
                bool dir = p2ID > p1ID;
                if (!dir) {
                  int tmp = p1ID;
                  p1ID = p2ID;
                  p2ID = tmp;
                }

                poreNE *p1 = poreIs[p1ID];
                poreNE *p2 = poreIs[p2ID];

                if (!p1)
                  std::cout << "  ERROR p1ID" << p1ID << std::endl;
                int tIDNext = throatIs.size();
                std::pair<std::map<int, int>::iterator, bool> ret =
                    p2->contacts.insert(std::pair<int, int>(p1ID, tIDNext));
                if (ret.second) {
                  if (ret.second !=
                      p1->contacts.insert(std::pair<int, int>(p2ID, tIDNext))
                          .second)
                    std::cout << "Errorpnm821-2 " << p1ID << " " << p2ID
                              << std::endl;
                  throatIs.push_back(new throatNE(tIDNext, p1ID, p2ID));
                  // if (p1ID<2 &&p2ID<2) cout <<"Errskziw2: p1ID<2 &&p2ID<2:
                  // tid: "<< tIDNext<<endl;
                }
                throatNE *trot = throatIs[ret.first->second];
                trot->CrosArea.y += (2 * dir - 1);
                // trot->C[1] += int3{{i-1,j-1,k-1}};
              }

              if (p1ID >= firstPore)
                ++(poreIs[p1ID]->surfaceArea);
              if (p2ID >= firstPore)
                ++(poreIs[p2ID]->surfaceArea);
            }
          }
        }

    (std::cout << throatIs.size() << " .. ").flush();

    for (int z = 0; z < VElems.nz() - 1; ++z) //     z >-<
      for (int y = 1; y < VElems.ny() - 1; ++y)
        for (int x = 1; x < VElems.nx() - 1; ++x) {
          int p1ID = VElems(x, y, z);
          if (p1ID >= 0) {
            int p2ID = VElems(x, y, z + 1);
            if (p1ID != p2ID) {
              if (p2ID >= 0) {
                bool dir = p2ID > p1ID;
                if (!dir) {
                  int tmp = p1ID;
                  p1ID = p2ID;
                  p2ID = tmp;
                }
                poreNE *p1 = poreIs[p1ID];
                poreNE *p2 = poreIs[p2ID];

                if (!p1)
                  std::cout << "  ERROR p1ID" << p1ID << std::endl;
                int tIDNext = throatIs.size();
                std::pair<std::map<int, int>::iterator, bool> ret =
                    p2->contacts.insert(std::pair<int, int>(p1ID, tIDNext));
                if (ret.second) {
                  if (ret.second !=
                      p1->contacts.insert(std::pair<int, int>(p2ID, tIDNext))
                          .second)
                    std::cout << "Errorpnm821-2 " << p1ID << " " << p2ID
                              << std::endl;
                  throatIs.push_back(new throatNE(tIDNext, p1ID, p2ID));
                  // if (p1ID<2 &&p2ID<2) cout <<"Errskziw3: p1ID<2 &&p2ID<2:
                  // tid: "<< tIDNext<<endl;
                }
                throatNE *trot = throatIs[ret.first->second];
                trot->CrosArea.z += (2 * dir - 1);
                // trot->C[2] += int3{{i-1,j-1,k-1}};
              }

              if (p1ID >= firstPore)
                ++(poreIs[p1ID]->surfaceArea);
              if (p2ID >= firstPore)
                ++(poreIs[p2ID]->surfaceArea);
            }
          }
        }

    (std::cout << throatIs.size() << ", ").flush();
  }

  std::cout << " nThroats: " << throatIs.size() << std::endl;

  nNodes = poreIs.size();
  nTrots = throatIs.size();
  nElems = nNodes + nTrots;

  std::cout << "nElems:  " << nElems << " = " << poreIs.size() << " + "
            << throatIs.size() << std::endl;

  throadAdditBalls.reserve(
      throatIs.size() * 5); // to improve the efficiency when later generating
                            // and storing new maximal balls for throat surfaces

  std::cout << "\ncalc throat properties: ";
  std::cout.flush();

  {
    std::cout << " collecting face voxels,  ";
    std::cout.flush();

    for (auto tr : throatIs) {
      tr->toxels2.reserve(std::abs(tr->CrosArea[0]) +
                          std::abs(tr->CrosArea[1]) +
                          std::abs(tr->CrosArea[2]) + 1);
      tr->toxels1.reserve(std::abs(tr->CrosArea[0]) +
                          std::abs(tr->CrosArea[1]) +
                          std::abs(tr->CrosArea[2]) + 1);
    }

    int nMultiTouchErrors = 0;
    forAllkji_1(VElems) {
      int p1ID = VElems(i, j, k);
      if (p1ID >= firstPore) {
        std::set<int> neis;
        int neiPID;

        neiPID = VElems(i - 1, j, k);
        if ((p1ID != neiPID) && (neiPID >= 0))
          neis.insert(neiPID);
        neiPID = VElems(i + 1, j, k);
        if ((p1ID != neiPID) && (neiPID >= 0))
          neis.insert(neiPID);
        neiPID = VElems(i, j - 1, k);
        if ((p1ID != neiPID) && (neiPID >= 0))
          neis.insert(neiPID);
        neiPID = VElems(i, j + 1, k);
        if ((p1ID != neiPID) && (neiPID >= 0))
          neis.insert(neiPID);
        neiPID = VElems(i, j, k - 1);
        if ((p1ID != neiPID) && (neiPID >= 0))
          neis.insert(neiPID);
        neiPID = VElems(i, j, k + 1);
        if ((p1ID != neiPID) && (neiPID >= 0))
          neis.insert(neiPID);
        for (int nei : neis) {
          poreNE *p1 = poreIs[p1ID];
          throatNE *trot = throatIs[p1->contacts[nei]];
          if (p1ID > nei) {
            dbgAsrt(vxlMap(k - 1, j - 1, i - 1));
            trot->toxels2.push_back(vxlMap(k - 1, j - 1, i - 1));
          } else {
            dbgAsrt(vxlMap(k - 1, j - 1, i - 1));
            trot->toxels1.push_back(vxlMap(k - 1, j - 1, i - 1));
          }
        }
        if (neis.size() > 1)
          ++nMultiTouchErrors;
      }
    }
    if (nMultiTouchErrors > 0)
      std::cout << "\n   Warning: " << nMultiTouchErrors
                << " voxels being in touch to more than two pores" << std::endl;
  }

  for (auto tr : throatIs) {
    if (tr->toxels2.empty() || tr->toxels2.empty())
      (std::cout << "  ERROR1017, toxls size:" << tr->toxels2.size() << " "
                 << tr->toxels1.size() << "  ")
          .flush();

    sort(tr->toxels2.begin(), tr->toxels2.end(), metaballcomparer());
    sort(tr->toxels1.begin(), tr->toxels1.end(), metaballcomparer());
  }

  std::cout << " calculating throat radii";
  std::cout.flush();
  for (throatNE *tr : throatIs) {
    if (tr->toxels2.size() > 0) {
      voxel *tvox2 =
          *(tr->toxels2.begin()); // get the largest distance map throat voxel
      if (tvox2->ball != nullptr) {
        medialBall *vbi = (*tr->toxels2.begin())->ball;
        vbi->type = 5;

        medialBall *mvi = vbi->mastrSphere();
        if (mvi != nullptr && mvi != vbi &&
            tr->e2 != VElems(mvi->fi + 1, mvi->fj + 1, mvi->fk + 1))
          std::cout << " Dmb2rrr  "
                    << VElems(mvi->fi + 1, mvi->fj + 1, mvi->fk + 1) << "   ";
      } else {
        tvox2->ball = new medialBall(tvox2, 15);
        throadAdditBalls.push_back(tvox2->ball);
        srf->moveUphill(tvox2->ball);
      }
    }

    if (!tr->toxels1.empty()) {
      sort(tr->toxels1.begin(), tr->toxels1.end(), metaballcomparer());
      voxel *tvox1 = *(tr->toxels1.begin());
      if (tvox1->ball) {
        medialBall *vbi = tvox1->ball;
        vbi->type = 6;

        medialBall *mvi = vbi->mastrSphere();
        if (mvi && mvi != vbi &&
            tr->e1 != VElems(mvi->fi + 1, mvi->fj + 1, mvi->fk + 1))
          std::cout << " Dmb1rrr  "
                    << VElems(mvi->fi + 1, mvi->fj + 1, mvi->fk + 1) << "   ";
      } else {
        tvox1->ball = new medialBall(tvox1, 16);
        throadAdditBalls.push_back(tvox1->ball);
        srf->moveUphill(tvox1->ball);
      }
    }
  }

  std::cout << "." << std::endl;
}

// #include "blockNet_vxlManip.cpp"
#include "blockNet_write_cnm.cpp"
