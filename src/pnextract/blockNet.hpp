#pragma once
#include "ElementGNE.hpp"
#include "medialSurf.hpp"
#include "poregrower.hpp"
#include "types.hpp"

#define _SideRadius 1

class blockNetwork {
public:
  blockNetwork(const ConfigSettings &cfg, BitVector3D &bool_img)
      : srf(cfg, bool_img), cfg(cfg), maxNCors(5), nBP6(cfg.nBP6) {};

  void createchildHierarchy();
  void CreateVElem();
  void createNewThroats();

public:
  medialSurface srf;
  // std::vector<mediaAxes*> mas;

  int firstPore;
  int lastPore;
  const ConfigSettings &cfg;

  std::array<voxel, 8> sides;

  TensorXXXDi32 VElems;

  std::vector<poreNE> poreIs;
  std::vector<throatNE> throatIs;

  std::vector<medialBall> poreAdditBalls;
  std::vector<medialBall> throadAdditBalls;

  int nNodes;
  int nTrots;
  int nElems;

  const int maxNCors;
  const int nBP6;
};
