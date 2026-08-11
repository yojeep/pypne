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

inline double randomG(std::mt19937 &rng) { /// to randomly distribute the shape
                                           /// factors, in case of errors
  // Box-Muller for standard normal distribution N(0,1)
  // original code: y = 0.00625 * (x1 * w + 5.0)
  //       = 0.00625 * x1*w  +  0.03125
  // where x1*w ~ N(0,1), then:
  //   mean = 0.03125, stddev = 0.00625
  static std::normal_distribution<double> dist(0.03125, 0.00625);

  double y = dist(rng);

  return (y > 0.049) ? 0.0625 : std::max(y, 0.01);
}

inline double randomG() {
  thread_local std::mt19937 rng{42};
  return randomG(rng);
}

class blockNetwork;
std::tuple<std::string, std::string, std::string, std::string>
get_network(blockNetwork &bn);
