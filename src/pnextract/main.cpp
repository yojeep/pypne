

#include "pybind11/eigen/tensor.h" // IWYU pragma: keep // Do not delete !!! The returned tensor is Eigen::Tensor
#include "pypne.hpp"
#include "blockNet.hpp"
#include "medialSurf.hpp"
#include <Eigen/CXX11/Tensor>
#include <pasta/bit_vector/bit_vector.hpp>


using namespace pybind11::literals;

auto genextraction(const int &nz, const int &ny, const int &nx,
                   const float &resolution, const py::array_t<uint8_t> &arr,
                   py::dict &config_dict, const int &n_workers) {
  num_workers = n_workers;
  auto buf = arr.request();
  bool *data = static_cast<bool *>(buf.ptr);

  size_t nvoxels = nz * ny * nx;
  pasta::BitVector bv_(nvoxels);

  for (size_t i = 0; i < nvoxels; ++i) {
    bv_[i] = data[i];
  }

  BitVector3D bool_img(bv_, nz, ny, nx);
  bool_img.buildRankQuery();

  auto cfg = ConfigSettings();

  cfg.read_config(config_dict);
  cfg.resolution = resolution, cfg.nz = nz, cfg.ny = ny, cfg.nx = nx;

  blockNetwork mpn(cfg, bool_img);
  mpn.srf.createBallsAndHierarchy();
  mpn.CreateVElem();
  mpn.createNewThroats();

  auto [link1, link2, node1, node2] = get_network(mpn);
  if (cfg.write_Statoil) {
    writeToFile((cfg.output_path + "_link1.dat"), link1);
    writeToFile((cfg.output_path + "_link2.dat"), link2);
    writeToFile((cfg.output_path + "_node1.dat"), node1);
    writeToFile((cfg.output_path + "_node2.dat"), node2);
  }

  py::dict pn("link1"_a = link1, "link2"_a = link2, "node1"_a = node1,
              "node2"_a = node2);

  py::dict pyres("VElems"_a = mpn.VElems, "pn"_a = pn);
  GlobalThreadPool::shutdown();
  return pyres;
}

PYBIND11_MODULE(pypne_cpp, m) { m.def("pnextract", &genextraction); }
