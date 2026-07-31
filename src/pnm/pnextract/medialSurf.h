#pragma once
#include "Eigen/CXX11/Tensor"
#include "inputData.h"
#include "typses.h"
#include "voxelmap.h"
#include <array>
#include <cstdint>
#include <mdspan>
#include <span>
#include <stddef.h>
#include <vector>
#include <voxelImage.h>

class medialSurface {
public:
  // class node {
  // public:
  //   inline constexpr node() : i(-32768), j(-32768), k(-32768) {};
  //   inline node(const node &v) : i(v.i), j(v.j), k(v.k) {};
  //   inline node(const voxel &v) : i(v.i), j(v.j), k(v.k) {};
  //   inline node(int ii, int jj, int kk) : i(ii), j(jj), k(kk) {};
  //   inline void operator=(const node &v) {
  //     i = v.i;
  //     j = v.j;
  //     k = v.k;
  //   }
  //   short i, j, k;
  // };

  medialSurface(inputDataNE &cfg);
  void setDefaults(double avgRad);
  void paradox_pre_removeincludedballI();
  void paradoxremoveincludedballI();

  void calc_distmaps();
  void buildvoxelspace();
  void smoothRadius();

  void competeForParentNoMerge(medialBall *vi, medialBall *vjv);
  void competeForParent(medialBall *vi, medialBall *vj);
  void findBoss(medialBall *);
  void createBallsAndHierarchy();

  void moveUphill(medialBall *b_i);
  void moveUphillp1(medialBall *b_i);

  voxel *vxl(int i, int j, int k) {
    if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz)
      return nullptr;

    segments &s = segs_[k * ny + j];
    int p = s.fsi(i);
    if (p != -1) {
      return (0 == s.s[p].value) ? (s.s[p].segV + (i - s.s[p].start)) : nullptr;
    }
    return nullptr;
  }

  const segment &segg(int i, int j, int k) const {
    if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz)
      return invalidSeg;

    const segments &s = segs_[k * ny + j];
    int p = s.fsi(i);
    if (p != -1)
      return (s.s[p]);
    std::cout << "Error can not find segment at " << i << " " << j << " " << k
              << " nSegs: " << s.cnt << std::endl;
    return (s.s[s.cnt]);
  }

  const segment &nextSegg(int i, int j, int k) const {
    if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz)
      return invalidSeg;

    const segments &s = segs_[k * ny + j];
    int p = s.fsi(i);
    if (p != -1)
      return (s.s[p + 1]);

    std::cout << "Error can not find next segment at " << i << " " << j << " "
              << k << " nSegs: " << s.cnt << std::endl;
    return (s.s[s.cnt]);
  }

  bool isInside(int i, int j, int k) const {
    return (i >= 0 && j >= 0 && k >= 0 && i < nx && j < ny && k < nz);
  }

  bool isJInside(int j) const { return (j >= 0 && j < ny); }

  bool isInside(int i) const { return (0 <= i && i < nx); }

public:
  const inputDataNE &cg_;
  int nx, ny, nz;
  size_t nVxls{0};
  size_t nBalls{0};

  std::vector<segments> &segs_;
  segment invalidSeg;
  std::vector<voxel> vxlSpace;
  Eigen::Tensor<voxel*, 3, Eigen::RowMajor> vxlMap;
  Eigen::TensorMap<Eigen::Tensor<uint8_t, 3, Eigen::RowMajor>>
      binary_image;
  Eigen::TensorMap<Eigen::Tensor<uint8_t, 1, Eigen::RowMajor>>
      binary_image_flat;

  std::vector<medialBall> ballSpace;
  medialBall ToBeAssigned{0};
  // std::vector<std::vector<size_t>> iZ;
  gtl::flat_hash_map<std::array<int, 3>, voxel> voxelmap;

  float _minRp;
  double _clipROutx;
  double _clipROutyz;
  double _midRf;
  double _MSNoise;
  double _lenNf;

  double _vmvRadRelNf;

  int _nRSmoothing;
  double _RCorsnf;
  float _RCorsn;
};