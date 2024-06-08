#ifndef __SYMMETRY_ANALYZER_PLANE_HPP__
#define __SYMMETRY_ANALYZER_PLANE_HPP__

#include <array>
#include <complex>
#include <set>
#include <tuple>
#include <variant>

#ifndef __EIGEN__
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
template <typename T>
using vector = std::vector<T, Eigen::aligned_allocator<T>>;
using vector3i = Eigen::Matrix<int, 3, 1>;
using vector3d = Eigen::Matrix<double, 3, 1>;
using matrix3d = Eigen::Matrix<double, 3, 3>;
using vector3cd = Eigen::Matrix<std::complex<double>, 3, 1>;
#endif

bool has_reduction(const vector3i& hkl);
inline bool has_reduction(int h, int k, int l) {
  return has_reduction(vector3i{h, k, l});
}

std::complex<double> structure_factor(const vector3d& pos, const vector3d& hkl);

std::tuple<vector3i, int> reduced_hkl(const vector3i& hkl);

class Planes {
 private:
  double d_min = 1.5;
  int extent = 6;
  double cp_factor;
  std::set<std::tuple<int, int, int>> planes;

 public:
  Planes(int extent, double d_min, double cp_factor);
  double get_cp_factor(const vector3i& hkl);
  vector<vector3d> search_close_packing_planes(int N_max = 10);

  double get_structure_factor(const vector3d hkl);
};

#endif
