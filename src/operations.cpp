#include "operations.hpp"

std::vector<vector3d> create_matrix(const PBC& pbc, const bool omit) {
  std::vector<int> i_list;
  std::vector<int> j_list;
  std::vector<int> k_list;

  std::vector<vector3d> retvector;

  if (pbc[0])
    i_list.emplace_back(-1, 0, 1);
  else
    i_list.emplace_back(0);

  if (pbc[1])
    j_list.emplace_back(-1, 0, 1);
  else
    j_list.emplace_back(0);

  if (pbc[2])
    k_list.emplace_back(-1, 0, 1);
  else
    k_list.emplace_back(0);

  for (const int i : i_list) {
    for (const int j : j_list) {
      for (const int k : k_list) {
	if (omit) {
	  if (i != 0 || j != 0 || k != 0)
	    retvector.push_back(vector3d{i, j, k});
	} else {
	  retvector.push_back(vector3d{i, j, k});
	}
      }
    }
  }
  return retvector;
}

std::vector<vector3d> filtered_coords(const std::vector<vector3d>& coords,
				      const PBC& pbc) {
  std::vector<vector3d> retvector = coords;
  for (size_t i = 0; i < coords.size(); i++) retvector[i] -= coords[i] * pbc[i];
  return retvector;
}
