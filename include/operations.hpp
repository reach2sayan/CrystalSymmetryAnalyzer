#ifndef __OPERATIONS_HPP__
#define __OPERATIONS_HPP__

#if !defined(__SYMMETRY_ANALYZER_OPERATIONS_HPP__)
#include "symmetry_operations.hpp"
#endif

#include <array>

using PBC = Eigen::Matrix<bool, 3, 1>;
using vector3d = Eigen::Matrix<double, 3, 1>;
using matrix3d = Eigen::Matrix<double, 3, 3>;

std::vector<vector3d> create_matrix(const PBC& pbc = {true, true, true},
				    bool omit = false);

std::vector<vector3d> filtered_coords(const std::vector<vector3d>& coords,
				      const PBC& pbc = {true, true, true});

#endif
