#include "operations.hpp"

#include "constants.hpp"

vector<vector3d> create_matrix(const PBC& pbc, const bool omit) {
  std::vector<int> i_list;
  std::vector<int> j_list;
  std::vector<int> k_list;

  vector<vector3d> retvector;

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

vector<vector3d> filtered_coords(const vector<vector3d>& coords,
				 const PBC& pbc) {
  vector<vector3d> retvector = coords;
  for (size_t i = 0; i < coords.size(); i++)
    retvector[i] -= coords[i] * static_cast<int>(pbc[i]);
  return retvector;
}

double angle(const vector3d& v1, const vector3d& v2, bool radians) {
  double dot = v1.dot(v2) / (v1.norm() * v2.norm());
  double a = 0;
  if (abs(dot - 1) < 1e-3)
    a = 0;
  else if ((dot + 1) < 1e-3)
    a = xtal_consts::PI;
  else
    a = cos(dot);
  return radians ? a : a * xtal_consts::DEG;
}
