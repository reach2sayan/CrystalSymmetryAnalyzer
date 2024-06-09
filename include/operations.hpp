#ifndef __OPERATIONS_HPP__
#define __OPERATIONS_HPP__

#if !defined(__SYMMETRY_ANALYZER_OPERATIONS_HPP__)
#include "symmetry_operations.hpp"
#endif

#include <array>

vector<vector3d> create_matrix(const PBC& pbc = {true, true, true},
			       bool omit = false);

vector<vector3d> filtered_coords(const vector<vector3d>& coords,
				 const PBC& pbc = {true, true, true});

#endif
