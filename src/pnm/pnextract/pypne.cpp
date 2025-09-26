#include <array>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "blockNet.h"
#include "writers.h"
#include "profilers.h"
#include <iostream>
#include "pypne.h"

namespace py = pybind11;
using namespace pybind11::literals;
using vec_str_pair = std::vector<std::pair<std::string, std::string>>;

vec_str_pair convert_dict(py::dict py_dict)
{
    vec_str_pair result;
    for (auto item : py_dict)
    {
        result.emplace_back(
            item.first.cast<std::string>(),
            item.second.cast<std::string>());
    }
    return result;
}

// template <typename T>
// py::array_t<T> vector_to_numpy(const std::vector<T> &vec)
// {
//     py::array_t<T> arr(vec.size());
//     auto buf = arr.request();
//     T *ptr = static_cast<T *>(buf.ptr);
//     std::copy(vec.begin(), vec.end(), ptr);
//     return arr;
// }

template <typename T>
py::array_t<T> vector_to_numpy(const std::vector<T> &vec)
{
    return py::array_t<T>(
        {vec.size()}, // shape
        {sizeof(T)},  // stride
        vec.data()    // 直接使用 vector 的内存
    );
}

auto genextraction(int nx, int ny, int nz, double resolution, py::array_t<unsigned char> arr, py::dict config_dict, int n_workers)
{
    num_workers = n_workers;
    auto buf = arr.request();
    std::vector<unsigned char> data(arr.size());
    // auto arr_uncheck = arr.unchecked<1>();
    // for (py::ssize_t i = 0; i < arr.size(); i++)
    // {
    //     data[i] = arr_uncheck(i);
    // }
    std::memcpy(data.data(), buf.ptr, arr.size() * sizeof(unsigned char));
    voxelImageT_PYPNE<unsigned char> vm;
    vm.setData(data);
    vm.setNnn(int3(nx, ny, nz));
    vm.setNij(nx * ny);
    std::string output_path = "pn";
    if (config_dict.contains("output_path"))
    {
        output_path = config_dict["output_path"].cast<string>();
    };

    auto cfg = inputDataNE_PYPNE(nx, ny, nz, resolution, vm, output_path);
    vec_str_pair cfg_data = convert_dict(config_dict);
    cfg.setData(cfg_data);
    cfg.createSegments();
    medialSurface *srf;
    blockNetwork mpn(srf, cfg);
    mpn.createMedialSurface(srf, cfg, 0);
    mpn.CreateVElem(0);

    mpn.createNewThroats(srf);
    if (cfg.getOr("write_Statoil", false))
        mpn.writePNM();

    if (cfg.getOr("write_hierarchy", false))
        vtuWriteMbMbs(cfg.name() + "_mbHierarchy" + _s(0), srf->ballSpace, mpn.poreIs, mpn.VElems, cfg.vxlSize, mpn.VElems.X0() + mpn.VElems.dx());
    if (cfg.getOr("write_throatHierarchy", false))
        vtuWriteThroatMbMbs(cfg.name() + "_throatHierarchy", mpn.throatIs, mpn.poreIs, mpn.VElems, cfg.vxlSize, mpn.VElems.X0() + mpn.VElems.dx());
    if (cfg.getOr("write_vtkNetwork", false))
        vtuWritePores(cfg.name() + "_pores", mpn.poreIs, mpn.throatIs, cfg.vxlSize, mpn.VElems.X0() + mpn.VElems.dx());
    if (cfg.getOr("write_vtkNetwork", false))
        vtuWriteTHroatSpheres(cfg.name() + "_throatsBalls", mpn.poreIs, mpn.throatIs, cfg.vxlSize, mpn.VElems.X0() + mpn.VElems.dx());
    int outputBlockSize = 0; /// keywords write_throats, write_poreMaxBalls and write_throatMaxBalls developed by Tom Bultreys
    if (cfg.giv("outputBlockSize", outputBlockSize))
        cout << "OutputBlockSize:" << outputBlockSize << endl;
    if (!outputBlockSize)
    {
        if (cfg.getOr("write_throats", false))
            VThroats(mpn).writeBin(cfg.name() + "_throats" + cfg.imgfrmt, 0, cfg.nx, 0, cfg.ny, 0, cfg.nz);
        if (cfg.getOr("write_poreMaxBalls", false))
            poreMaxBalls(mpn).writeBin(cfg.name() + "_poreMBs" + cfg.imgfrmt, 0, cfg.nx, 0, cfg.ny, 0, cfg.nz);
        if (cfg.getOr("write_throatMaxBalls", false))
            throatMaxBalls(mpn).writeBin(cfg.name() + "_throatMBs" + cfg.imgfrmt, 0, cfg.nx, 0, cfg.ny, 0, cfg.nz);
    }
    else
    {
        int blockNumber = 1, beginSlice = 0, endSlice = 0;
        while (endSlice < cfg.nz - 1)
        {
            cout << " WRITING BLOCK \n";
            beginSlice = (blockNumber - 1) * outputBlockSize;
            endSlice = min(blockNumber * outputBlockSize, cfg.nx - 1);
            if (cfg.getOr("write_throats", false))
                VThroats(mpn, beginSlice, endSlice).writeBin(cfg.name() + "_throats" + _s(blockNumber) + cfg.imgfrmt, 0, cfg.nx, 0, cfg.ny, 0, endSlice - beginSlice);
            if (cfg.getOr("write_poreMaxBalls", false))
                poreMaxBalls(mpn, beginSlice, endSlice).writeBin(cfg.name() + "_poreMBs" + _s(blockNumber) + cfg.imgfrmt, 0, cfg.nx, 0, cfg.ny, 0, endSlice - beginSlice);
            if (cfg.getOr("write_throatMaxBalls", false))
                throatMaxBalls(mpn, beginSlice, endSlice).writeBin(cfg.name() + "_throatMBs" + _s(blockNumber) + cfg.imgfrmt, 0, cfg.nx, 0, cfg.ny, 0, endSlice - beginSlice);
            blockNumber++;
        }
    }

    if (cfg.getOr("write_vtkNetwork", false))
        vtuWriteThroats(cfg.name() + "_throats", mpn.poreIs, mpn.throatIs, cfg.vxlSize, mpn.VElems.X0() + mpn.VElems.dx());

    cout << endl
         << cfg.name() << endl;
    cout << "***  " << mpn.poreIs.size() << "-" << mpn.nBP6 << " pores, " << mpn.throatIs.size() << " throats,   ratio: " << double(mpn.throatIs.size()) / (mpn.poreIs.size() - 6. + 1e-6) << "  ***" << endl;
    cout << "end" << endl;

    auto VElems = mpn.VElems.data_;
    py::array_t<int> py_VElems = vector_to_numpy<int>(VElems);

    auto [link1, link2, node1, node2] = compute_PYPNEData_string(mpn);
    py::dict pn("link1"_a = link1, "link2"_a = link2, "node1"_a = node1, "node2"_a = node2);
    py::dict pyres("VElems"_a = py_VElems, "pn"_a = pn);
    GlobalThreadPool::shutdown();
    return pyres;
}

PYBIND11_MODULE(pypne_cpp, m)
{
    m.def("pnextract", &genextraction);
}
