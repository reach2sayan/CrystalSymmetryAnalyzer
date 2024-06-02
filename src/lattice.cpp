#include "lattice.hpp"

Lattice::Lattice(const LType _ltype, const double _volume,
		 const matrix3d& _matrix, const PBC& _pbc, bool _random,
		 bool _allow_volume_reset, const LUniqAx _uax)
    : pbc{_pbc},
      ltype{_ltype},
      random{_random},
      allow_volume_reset{_allow_volume_reset},
      unique_axis{_uax},
      matrix{_matrix},
      volume{_volume} {
  dim = std::accumulate(std::begin(pbc), std::end(pbc), 0);

  switch (_ltype) {
    case LType::triclinic:
      norm_matrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
      break;
    case LType::monoclinic:
      if (_pbc == std::array<int, 3>{1, 1, 1}) {
	norm_matrix << 1, 0, 0, 0, 1, 0, 1, 0, 1;
      } else {
	switch (_uax) {
	  case LUniqAx::a:
	    norm_matrix << 1, 0, 0, 0, 1, 0, 0, 1, 1;
	    break;
	  case LUniqAx::b:
	    norm_matrix << 1, 0, 0, 0, 1, 0, 1, 0, 1;
	    break;
	  case LUniqAx::c:
	    norm_matrix << 1, 0, 0, 1, 1, 0, 0, 0, 1;
	    break;
	}
      }
      break;

    case LType::orthorhombic:
    case LType::tetragonal:
    case LType::trigonal:
    case LType::hexagonal:
    case LType::cubic:
      norm_matrix = Eigen::Matrix3d::Identity();
      break;
    case LType::spherical:
    case LType::ellipsoidal:
      norm_matrix = Eigen::Matrix3d::Zero();
      break;
  }
  stress_normalization_matrix = norm_matrix;

  switch (_ltype) {
    case LType::tetragonal:
    case LType::trigonal:
    case LType::hexagonal:
      stress_indices.reserve(2);
      stress_indices.emplace_back(std::pair<int, int>(0, 0));
      stress_indices.emplace_back(std::pair<int, int>(1, 1));
      break;
    case LType::cubic:
      stress_indices.reserve(3);
      stress_indices.emplace_back(std::pair<int, int>(0, 0));
      stress_indices.emplace_back(std::pair<int, int>(1, 1));
      stress_indices.emplace_back(std::pair<int, int>(2, 2));
      break;
    default:
      stress_indices;
      break;
  }

  if (_ltype == LType::triclinic)
    a_tol = 15.0;
  else
    a_tol = 9.9;
  dim = std::accumulate(std::begin(_pbc), std::end(_pbc), 0);
  dof = get_dofs();
}

int Lattice::get_dofs() const {
  switch (ltype) {
    case LType::triclinic:
      return 6;
    case LType::monoclinic:
      return 4;
    case LType::orthorhombic:
      return 3;
    case LType::tetragonal:
    case LType::hexagonal:
    case LType::trigonal:
      return 2;
    default:
      return 1;
  }
}

