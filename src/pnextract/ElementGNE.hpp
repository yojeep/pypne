#pragma once
#include "globals.hpp"
#include "gtl/phmap.hpp"
#include "types.hpp"

inline int nAprox(1);
class medialBall;

class voxel {
public:
  voxel(int z, int y, int x, float rR)
      : z(z), y(y), x(x), R(rR), ball(nullptr) {};
  voxel() : z(-1), y(-1), x(-1), R(0.f), ball(nullptr) {};
  explicit operator bool() const noexcept { return z != -1; }

public:
  int z, y, x;
  float R;
  medialBall *ball;
};

/// maximal-sphere class, describing spheres, which are generated for voxels on
/// the medial surface
class medialBall {
public:
  float fz, fy, fx;
  float R;
  medialBall *boss;
  voxel *vxl;
  unsigned short corId;
  short type;

public:
  medialBall() = delete;
  medialBall(short t)
      : fz(-10000.0f), fy(-0.5f), fx(-10000.0f), R(-10000.0f), boss(this),
        vxl(nullptr), corId(0), type(t) {}
  medialBall(voxel *v, short t)
      : fz(v->z + _0p5), fy(v->y + _0p5), fx(v->x + _0p5), R(v->R), boss(this),
        vxl(v), corId(0), type(t) {}
  ~medialBall() {}

  int iz() const { return static_cast<int>(fz); }
  int iy() const { return static_cast<int>(fy); }
  int ix() const { return static_cast<int>(fx); }

  const int level() const {
    const medialBall *current = this;
    int count = 1;
    while (current != current->boss) {
      medialBall *nxt = current->boss;
      __builtin_prefetch(nxt->boss, 0, 1);
      current = nxt;
      ++count;
    }
    return count;
  }

  medialBall *mastrSphere() {
    medialBall *current = this;
    while (current != current->boss) {
      medialBall *nxt = current->boss;
      __builtin_prefetch(nxt->boss, 0, 1);
      current = nxt;
    }
    return current;
  }

  const medialBall *mastrSphere() const {
    const medialBall *current = this;
    while (current != current->boss) {
      medialBall *nxt = current->boss;
      __builtin_prefetch(nxt->boss, 0, 1);
      current = nxt;
    }
    return current;
  }

  const bool inParents(const medialBall *vj) const {

    const medialBall *current = boss;
    while (current != current->boss) {
      medialBall *nxt = current->boss;
      __builtin_prefetch(nxt->boss, 0, 1);
      if (current == vj) {
        return true;
      }
      current = nxt;
    }
    return current == vj;
  }
};

class poreNE {

public:
  // poreNE(const poreNE &) {};
  poreNE(const poreNE &) = delete;
  poreNE &operator=(const poreNE &) = delete;
  poreNE(poreNE &&other) noexcept = default;
  poreNE &operator=(poreNE &&other) noexcept = default;

  poreNE() : surfaceArea(0), volumn(0), mb(nullptr) {}
  poreNE(int sa, int volumn, const medialBall *mb)
      : surfaceArea(sa), volumn(volumn), mb(mb) {}
  // dbl3 node() const { return mb->node(); }
  float radius() const { return mb->R; }

public:
  int surfaceArea;
  int volumn;

  gtl::flat_hash_map<int, int> contacts;
  const medialBall *mb;
};

class throatNE {
public:
  throatNE(const throatNE &) = delete;
  throatNE &operator=(const throatNE &) = delete;
  throatNE(throatNE &&other) noexcept = default;
  throatNE &operator=(throatNE &&other) noexcept = default;

  throatNE(int tid, int elm1, int elm2)
      : tid(tid), e1(elm1), e2(elm2), surfaceArea(0), volumn(0), cachInd(0),
        nCrnrs(0), CrosArea{0., 0., 0.} {}
  int nToxel2Balls() const {
    int nBs = 0;
    for (auto vi : toxels2) {
      if (vi->ball)
        ++nBs;
    }
    return nBs;
  }

  float radius() const {
    return toxels1.empty() ? mb22()->R : 0.5 * (mb22()->R + mb11()->R);
  }

public:
  int tid;
  int e1;
  int e2;

  int surfaceArea;
  int volumn;

  unsigned short cachInd;

  short nCrnrs;

  // std::vector<vxlface>  vxfaces;
  std::vector<voxel *> toxels1;
  std::vector<voxel *> toxels2;

  const medialBall *mb11() const {
    return toxels1.empty() ? nullptr : toxels1.front()->ball;
  }
  const medialBall *mb22() const { return toxels2.front()->ball; }

  // dbl3 CrosArea;
  Vector3f32 CrosArea;
  // int3x3 C;
};

inline float distSq(const medialBall *i, const medialBall *j) {
  const float dz = i->fz - j->fz;
  const float dy = i->fy - j->fy;
  const float dx = i->fx - j->fx;
  return dz * dz + dy * dy + dx * dx;
}

inline float distSq(const medialBall &i, const medialBall &j) {
  const float dz = i.fz - j.fz;
  const float dy = i.fy - j.fy;
  const float dx = i.fx - j.fx;
  return dz * dz + dy * dy + dx * dx;
}

inline float distSq(const voxel &i, const voxel &j) {
  const int dx = i.x - j.x;
  const int dy = i.y - j.y;
  const int dz = i.z - j.z;
  return dx * dx + dy * dy + dz * dz;
}

class metaballcomparer {
public:
  bool operator()(const voxel *a, const voxel *b) const {
    return (a->R) > (b->R);
  }
};

inline const medialBall &pos_of(const medialBall &b) { return b; }
inline const medialBall &pos_of(const medialBall *b) { return *b; }

template <typename A, typename B> inline float dist(const A &a, const B &b) {
  const auto &aa = pos_of(a);
  const auto &bb = pos_of(b);

  const float dx = aa.fx - bb.fx;
  const float dy = aa.fy - bb.fy;
  const float dz = aa.fz - bb.fz;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}