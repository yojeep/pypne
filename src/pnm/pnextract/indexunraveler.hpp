#pragma once
#include <vector>

class IndexUnraveler {
private:
  std::vector<size_t> strides_;
  std::vector<size_t> shape_;

public:
  explicit IndexUnraveler(const std::vector<size_t> &shape)
      : shape_(shape.begin(), shape.end()) {
    const size_t ndim = shape_.size();
    strides_.resize(ndim);
    if (ndim > 0) {
      strides_.back() = 1;
      for (size_t i = ndim - 1; i > 0; --i) {
        strides_[i - 1] = strides_[i] * shape_[i];
      }
    }
  }

  template <typename CoordsT>
  void unravel(size_t flat_idx, CoordsT &coords) const {
    // assert(coords.size() == shape_.size() && "Coords size must match shape");

    size_t remainder = flat_idx;
    const size_t n = shape_.size();
    for (size_t i = 0; i < n; ++i) {
      coords[i] = remainder / strides_[i];
      remainder %= strides_[i];
    }
  }
};