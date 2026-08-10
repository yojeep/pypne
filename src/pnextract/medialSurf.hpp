#pragma once
#include <vector>
#include <cstddef>
#include "config.hpp"
#include "ElementGNE.hpp"
#include "Eigen/CXX11/Tensor"
#include "bit_vector.hpp"

class medialSurface {
public:
  medialSurface(const ConfigSettings &cfg, BitVector3D &bool_img)
      : cfg(cfg), nz(cfg.nz), ny(cfg.ny), nx(cfg.nx), bool_img(bool_img),
        vxlMap(bool_img.bv_, nz, ny, nx) {
    nVxls = bool_img.count1();
    if (!nVxls) {
      std::cout << " no voxels no balls,\n" << std::endl;
      exit(1);
    }
    nvoxels = bool_img.size();
    vxlMap.buildRankQuery();
    setDefaults(-5.);
  }
  void buildvoxelspace();
  float calc_distmaps();
  void setDefaults(float avgRad);
  void paradox_pre_removeincludedballI();
  void paradoxremoveincludedballI();
  void smoothRadius();
  void competeForParent(medialBall *vi, medialBall *vj);
  void findBoss(medialBall *vi);
  void createBallsAndHierarchy();
  void moveUphill(medialBall &bi);
  void moveUphillp1(medialBall &bi);

public:
  const ConfigSettings &cfg;
  const int nz, ny, nx;
  const BitVector3D &bool_img;
  size_t nVxls{0};
  size_t nvoxels{0};
  size_t nBalls{0};

  BitMap3D<voxel> vxlMap;

  std::vector<medialBall> ballSpace;
  medialBall ToBeAssigned{0};

  float _minRp;
  float _clipROutx;
  float _clipROutyz;
  float _midRf;
  float _MSNoise;
  float _lenNf;
  float _vmvRadRelNf;
  int _nRSmoothing;
  float _RCorsnf;
  float _RCorsn;
};