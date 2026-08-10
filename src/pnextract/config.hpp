#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <optional>
#include <string>

namespace py = pybind11;
class ConfigSettings {

private:
  template <typename T>
  void assign_if_exists(const py::dict &cfg, const char *key, T &target) {
    if (!cfg[key].is_none()) {
      target = cfg[key].cast<T>();
    }
  }
#define READ_CFG(key) assign_if_exists(cfg, #key, key)

public:
  float resolution;
  int nx;
  int ny;
  int nz;
  int nBP6 = 2;
  std::optional<float> minRp;
  std::optional<float> clipROutx;
  std::optional<float> clipROutyz;
  std::optional<float> midRf;
  std::optional<float> MSNoise;
  std::optional<float> lenNf;
  std::optional<float> vmvRadRelNf;
  std::optional<int> nRSmoothing;
  std::optional<float> RCorsnf;
  std::optional<float> RCorsn;
  bool write_Statoil = false;
  bool write_elements = false;
  bool write_all = false;
  std::string output_path = "";

  void read_config(py::dict &cfg) {
    READ_CFG(nBP6);
    READ_CFG(minRp);
    READ_CFG(clipROutx);
    READ_CFG(clipROutyz);
    READ_CFG(midRf);
    READ_CFG(MSNoise);
    READ_CFG(lenNf);
    READ_CFG(vmvRadRelNf);
    READ_CFG(nRSmoothing);
    READ_CFG(RCorsnf);
    READ_CFG(RCorsn);
    READ_CFG(write_Statoil);
    READ_CFG(write_elements);
    READ_CFG(write_all);
    READ_CFG(output_path);
  };
};
