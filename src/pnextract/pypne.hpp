#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "blockNet.hpp"
#include <fstream>
namespace py = pybind11;

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
py::array_t<T> vector_to_numpy(const std::vector<T> &vec) {
  return py::array_t<T>({vec.size()}, // shape
                        {sizeof(T)},  // stride
                        vec.data()    // 直接使用 vector 的内存
  );
}

inline bool writeToFile(const std::string &filename,
                        const std::string &content) {
  std::ofstream outFile(filename, std::ios::out);
  if (outFile.is_open()) {
    outFile << content;
    outFile.close();
    return true;
  }
  std::cerr << "Error: Cannot open file for writing: " << filename << std::endl;
  return false;
}

// #include <format>

inline double
randomG() /// to randomly distribute the shape factors, in case of errors
{
  double x1, x2, w, y;
  do {
    x1 = 2. * rand() / RAND_MAX - 1.;
    x2 = 2. * rand() / RAND_MAX - 1.;
    w = x1 * x1 + x2 * x2;
  } while (w >= 1.);
  w = sqrt((-2. * log(w)) / w);
  y = 0.00625 * (x1 * w + 5.);
  if (y > 0.049)
    y = 0.0625;
  return y;
}

class blockNetwork;
std::tuple<std::string, std::string, std::string, std::string>
get_network(blockNetwork &bn);
