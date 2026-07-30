#include "blockNet.h"
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

  VElems.reset(cg.nx + 2, cg.ny + 2, cg.nz + 2, -257);
  VElems.X0Ch() = cg.VImage.X0() - cg.VImage.dx();
  VElems.dxCh() = cg.VImage.dx();
  int nVVs = 0;
  for (int iz = 0; iz < cg.nz; ++iz) {
    for (int iy = 0; iy < cg.ny; ++iy) {
      const segments &s = cg.segs_[iz * cg.ny + iy];
      for (int ix = 0; ix < s.cnt; ++ix) {
        int value = -1 - int(s.s[ix].value);
        nVVs = std::max(nVVs, value);
        for (int i = s.s[ix].start; i < s.s[ix + 1].start; ++i) {
          VElems(i + 1, iy + 1, iz + 1) = value;
        }
      }
    }
  }
  ++nVVs;

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
                                    // be sued to uniquely identify pores
    VElems.setSlice('j', cg.ny + 1, -1 - 3 - nVVs);
    VElems.setSlice('k', 0, -1 - 4 - nVVs);
    VElems.setSlice('k', cg.nz + 1, -1 - 5 - nVVs);
  }

  firstPore = 2;

  {

    int uasyned = -1;

    const std::vector<medialBall> &balspc = srf->ballSpace;

    firstPores = (poreIs.size());
    for (const auto &bi : balspc) {
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
    auto &binary_image = cg.VImage.data_;
    for (const auto &bi : balspc) {
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
                z >= cg.nz ||
                binary_image[z * cg.ny * cg.nx + y * cg.nx + x] != 0)
              continue;
            int zVE = z + 1;
            int yVE = y + 1;
            int xVE = x + 1;
            int Vi = VElems(xVE, yVE, zVE);
            if (Vi == uasyned)
              VElems(xVE, yVE, zVE) = Vm;
            else if (Vi != Vm) {
              voxel *vi = srf->vxl(x, y, z);
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
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    std::cout << std::endl;
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    std::cout << std::endl;
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedStrict(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);
    std::cout << std::endl;
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    std::cout << std::endl;
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);

    std::cout << std::endl;

    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    std::cout << std::endl;
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    std::cout << std::endl;
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    std::cout << std::endl;
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores_X2(VElems, voxls, firstPores, lastPores, uasyned);
    growPores_X2(VElems, voxls, firstPores, lastPores, uasyned);
    growPores_X2(VElems, voxls, firstPores, lastPores, uasyned);
    growPores_X2(VElems, voxls, firstPores, lastPores, uasyned);
    std::cout << std::endl;

    medianElem(cg, VElems, voxls, firstPores, lastPores, poreIs);
    medianElem(cg, VElems, voxls, firstPores, lastPores, poreIs);
    medianElem(cg, VElems, voxls, firstPores, lastPores, poreIs);

    medianElem(cg, VElems, voxls, firstPores, lastPores, poreIs);
    medianElem(cg, VElems, voxls, firstPores, lastPores, poreIs);

    std::cout << std::endl;

    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    while (growPores_X2(VElems, voxls, firstPores, lastPores, uasyned))
      ;
    growPores(VElems, voxls, firstPores, lastPores, uasyned);

    std::cout << std::endl;

    retreatPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs,
                       uasyned);

    for (const auto &bi : balspc) {
      medialBall *mastrSphere = bi.mastrSphere();
      VElems(bi.fi + 1, bi.fj + 1, bi.fk + 1) =
          VElems(mastrSphere->fi + 1, mastrSphere->fj + 1, mastrSphere->fk + 1);
    }
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedian(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores_X2(VElems, voxls, firstPores, lastPores, uasyned);
    growPoresMedEqs(cg, VElems, voxls, firstPores, lastPores, poreIs, uasyned);
    growPores(VElems, voxls, firstPores, lastPores, uasyned);
    growPores_X2(VElems, voxls, firstPores, lastPores, uasyned);

    std::cout << std::endl;

    medianElem(cg, VElems, voxls, firstPores, lastPores, poreIs);
    for (const auto &bi : balspc) {
      medialBall *mastrSphere = bi.mastrSphere();
      VElems(bi.fi + 1, bi.fj + 1, bi.fk + 1) =
          VElems(mastrSphere->fi + 1, mastrSphere->fj + 1, mastrSphere->fk + 1);
    }
    growPoresMedEqsLoose(cg, VElems, voxls, firstPores, lastPores, poreIs,
                         uasyned);
    int label = firstPore;
    for (const auto &bi : balspc) {
      if (bi.boss == &bi) {
        VElems(bi.fi + 1, bi.fj + 1, bi.fk + 1) = label;
        ++label;
      }
    }

    std::cout << std::endl;
  }
}

void blockNetwork::createNewThroats(medialSurface *&srf) {

  // vector<int> iThroatFaces;
  // iThroatFaces.reserve(throatIs.size()+poreIs.size()*10);
  {
    std::cout << "\nlooking for connections, iFirst = " << throatIs.size()
              << " ..  ";
    std::cout.flush();

    for (int i = 0; i < VElems.nx() - 1; i++)
      for (int k = 1; k < VElems.nz() - 1; ++k)
        for (int j = 1; j < VElems.ny() - 1; ++j) {
          int p1ID = VElems(i, j, k);
          if (p1ID >= 0) {
            int p2ID = VElems(i + 1, j, k);
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

    for (int j = 0; j < VElems.ny() - 1; j++) //     y >-<
      for (int k = 1; k < VElems.nz() - 1; ++k)
        for (int i = 1; i < VElems.nx() - 1; ++i) {
          int p1ID = VElems(i, j, k);
          if (p1ID >= 0) {
            int p2ID = VElems(i, j + 1, k);
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

    for (int k = 0; k < VElems.nz() - 1; ++k) //     z >-<
      for (int j = 1; j < VElems.ny() - 1; ++j)
        for (int i = 1; i < VElems.nx() - 1; ++i) {
          int p1ID = VElems(i, j, k);
          if (p1ID >= 0) {
            int p2ID = VElems(i, j, k + 1);
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
            dbgAsrt(srf->vxl(i - 1, j - 1, k - 1));
            trot->toxels2.push_back((srf->vxl(i - 1, j - 1, k - 1)));
          } else {
            dbgAsrt(srf->vxl(i - 1, j - 1, k - 1));
            trot->toxels1.push_back((srf->vxl(i - 1, j - 1, k - 1)));
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

#include "blockNet_vxlManip.cpp"
#include "blockNet_write_cnm.cpp"
