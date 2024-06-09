#include "plane.hpp"

#include <algorithm>
#include <numeric>

#ifndef __SYMMETRY_ANALYZER_CONSTANTS_HPP__
#include "constants.hpp"
#endif

bool has_reduction(const vector3i& hkl) {
  int gcf = std::__gcd(hkl[0], std::__gcd(hkl[1], hkl[2]));
  auto it = std::find_if_not(std::begin(hkl), std::end(hkl),
			     [](int m) { return m == 0; });
  auto first_zero_index = std::distance(std::begin(hkl), it);
  if (gcf > 1)
    return true;
  else if (first_zero_index < 3 && hkl[first_zero_index] < 0)
    return true;
  else
    return false;
}

std::complex<double> structure_factor(const vector3d& pos,
				      const vector3d& hkl) {
  return std::exp(-2.0 * xtal_consts::PI * std::complex<double>{0, 1} *
		  pos.dot(hkl));
}

std::tuple<vector3i, int> reduced_hkl(const vector3i& hkl) {
  int gcf = std::__gcd(hkl[0], std::__gcd(hkl[1], hkl[2]));

  if (gcf > 1) {
    return {hkl / gcf, gcf};
  } else {
    return {hkl, 1};
  }
}

double d_spacing(const matrix3d& inverted_matrix, const vector3d& hkl) {
  return 1 / (inverted_matrix * hkl).norm();
}

