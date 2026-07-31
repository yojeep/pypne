#pragma once

#include "ElementGNE.h"
#include "inputData.h"
#include "medialSurf.h"
#include "poregrower.hpp"
#include <cstddef>
#include <cstdint>
#include <vector>

#define _SideRadius 1

class blockNetwork {
public:
  blockNetwork(medialSurface *&rs, const inputDataNE &cfg)
      : srf(rs), cg(cfg), maxNCors(5) {
    nBP6 = cfg.nBP6;
    SideImax = nBP6 - 1;
  };

  void createMedialSurface(medialSurface *&srf, inputDataNE &cfg,
                           size_t startValue);
  void collectAllballs(medialSurface *&srf, inputDataNE &cfg,
                       size_t startValue);
  void createchildHierarchy();
  void CreateVElem(size_t startValue);
  void createNewThroats(medialSurface *&srf);
  void findParentInElem(medialBall *vi, int e1s2, int &EbadElem, int &EbadBoss,
                        int &nAllElem, int &nAllBoss);
  void createThroatBallConnectivity(medialSurface *&srf);
  void createCornerHirarchy();

  void writePNM() const;

  static int SideImax;
  static int nBP6; // 1+SideImax

public:
  medialSurface *&srf;
  // std::vector<mediaAxes*> mas;
  int firstPore;
  int firstPores;
  int lastPores;
  const inputDataNE &cg;
  std::array<voxel, 8> sides;

  voxelImageT<int> VElems;

  std::vector<poreNE *> poreIs;
  std::vector<throatNE *> throatIs;

  // Eigen::TensorMap<Eigen::Tensor<const uint8_t, 3, Eigen::RowMajor>>
  //     &binary_image;
  // std::vector<voxel> vxlSpace;
  // const Eigen::Tensor<voxel *, 3, Eigen::RowMajor> vxlMap;
  // std::vector<medialBall> ballSpace;

  std::vector<medialBall *> throadAdditBalls;

  int nNodes;
  int nTrots;
  int nElems;

  const int maxNCors;
};
