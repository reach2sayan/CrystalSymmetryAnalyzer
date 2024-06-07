#include "spherical_harmonics.hpp"

#include <cassert>
#include <cmath>
#include <vector>
#if !defined(__SYMMETRY_ANALYZER_CONSTANTS_HPP__)
#include "constants.hpp"
#endif

#ifndef __EIGEN__
#include <eigen3/Eigen/Dense>
template <typename T>
using vector = std::vector<T, Eigen::aligned_allocator<T>>;
using matrix4d = Eigen::Matrix<double, 4, 4>;
using vector4d = Eigen::Matrix<double, 4, 1>;
using matrix3d = Eigen::Matrix<double, 3, 3>;
using vector3d = Eigen::Matrix<double, 3, 1>;
#endif

constexpr double factorial(int val) {
  if (val <= 16) {
    return FACTORIAL_CACHE[val];
  } else {
    return val * factorial(val - 1);
  }
}

double spherical_harmonics(int l, int m, double phi, double theta) {
  assert(l >= 0);
  assert(-l <= m && m <= l);

  double kml = sqrt((2.0 * l + 1) * factorial(l - abs(m)) /
		    (4.0 * xtal_consts::PI * factorial(l + abs(m))));
  if (m > 0) {
    return sqrt(2.0) * kml * cos(m * phi) * std::sph_legendre(l, m, cos(theta));
  } else if (m < 0) {
    return sqrt(2.0) * kml * sin(-m * phi) *
	   std::sph_legendre(l, -m, cos(theta));
  } else {
    return kml * std::sph_legendre(l, 0, cos(theta));
  }
}
