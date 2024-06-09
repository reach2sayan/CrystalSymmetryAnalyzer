#include "lattice.hpp"

#include <cmath>
#include <random>

#include "constants.hpp"

#if !defined(__OPERATIONS_HPP__)
#include "operations.hpp"
#endif

#if !defined(__SYMMETRY_ANALYZER_OPERATIONS_HPP__)
#include "symmetry_operations.hpp"
#endif

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
      if (_pbc[0] == 1 && _pbc[1] == 1 && _pbc[2] == 1) {
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
      stress_indices = {};
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

std::pair<vector<vector3d>, std::vector<double>> Lattice::get_lengths() const {
  vector<vector3d> mat = create_matrix(pbc);
  vector<vector3d> retvector;
  std::vector<double> retnorm;

  retvector.reserve(mat.size());
  retnorm.reserve(mat.size());

  for (const vector3d& v : mat) {
    retvector.emplace_back(matrix * v);
    retnorm.push_back((matrix * v).norm());
  }

  return {retvector, retnorm};
}

CellParams Lattice::get_parameters(bool degree) const {
  double _alpha = alpha;
  double _gamma = gamma;
  double _beta = beta;
  if (degree) {
    _alpha = xtal_consts::DEG * alpha;
    _beta = xtal_consts::DEG * beta;
    _gamma = xtal_consts::DEG * gamma;
  };

  return {a, b, c, _alpha, _beta, _gamma};
}

double Lattice::get_worst_angle() const {
  return std::max({fabs(alpha - xtal_consts::PI / 2),
		   fabs(gamma - xtal_consts::PI / 2),
		   fabs(beta - xtal_consts::PI / 2)});
}

vector<matrix3d> Lattice::get_permutation_matrices() const {
  vector<matrix3d> retvector;
  switch (ltype) {
    case LatticeType::monoclinic:
      retvector.reserve(2);
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 0, 0, 1, 0, 1, 0, 1, 0, 0).finished());
      break;
    case LatticeType::triclinic:
      retvector.reserve(4);
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 0, 1, 0, 1, 0).finished());
      retvector.emplace_back(
	  (matrix3d() << 0, 0, 1, 0, 1, 0, 1, 0, 0).finished());
      retvector.emplace_back(
	  (matrix3d() << 0, 1, 0, 1, 0, 0, 0, 0, 1).finished());
    default:
      retvector.emplace_back(matrix3d::Identity());
  }
  return retvector;
}

vector<matrix3d> Lattice::get_transformation_matrices() const {
  vector<matrix3d> retvector;
  switch (ltype) {
    case LatticeType::monoclinic:
      retvector.reserve(6);
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, 1, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, -1, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 1, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, -1, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished());
      break;
    case LatticeType::triclinic:
      retvector.reserve(16);
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, -1, 0, 0, 0, -1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, 1, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, -1, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 1, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, -1, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, 0, 1, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 1, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, 0, -1, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, -1, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 1, 0, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, -1, 0, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 1, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, -1, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << -1, 0, 0, 0, 1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, -1, 0, 0, 0, 1).finished());
      retvector.emplace_back(
	  (matrix3d() << 1, 0, 0, 0, 1, 0, 0, 0, -1).finished());
      break;
    default:
      retvector.emplace_back(matrix3d::Identity());
      break;
  }
  return retvector;
}

CellParams generate_cell_parameters(const LatticeType ltype, const double vol,
				    const double minvec, const double minangle,
				    const double max_ratio,
				    const int max_attempts) {
  const double maxangle = xtal_consts::PI - minangle;
  int attempts = 0;
  while (attempts < max_attempts) {
    double a = 0;
    double b = 0;
    double c = 0;

    switch (ltype) {
      case LatticeType::triclinic:

	break;
    }
  }
  return {1, 1, 1, 90, 90, 90};
}

matrix3d random_shear_matrix(double width, bool unitary) {
  matrix3d retmatrix = matrix3d::Identity();
  double determinant = 0;
  std::normal_distribution<double> rnd{0, width};
  while (determinant == 0) {
    //[[1, a, b], [a, 1, c], [b, c, 1]])
    retmatrix(0, 1) = retmatrix(1, 0) = rnd(rand_generator);
    retmatrix(0, 2) = retmatrix(2, 0) = rnd(rand_generator);
    retmatrix(1, 2) = retmatrix(2, 1) = rnd(rand_generator);
    determinant = retmatrix.determinant();
  }
  if (unitary) retmatrix /= cbrt(determinant);
  return retmatrix;
}

vector3d random_vector(const vector3d& minvector, const vector3d& maxvector,
		       double width, bool unit) {
  std::normal_distribution<double> rnd{0, width};
  vector3d retvec = {exp(rnd(rand_generator)), exp(rnd(rand_generator)),
		     exp(rnd(rand_generator))};
  if (unit) retvec /= retvec.norm();
  return retvec;
}

/*
def random_vector(minvec=[0.0, 0.0, 0.0], maxvec=[1.0, 1.0, 1.0], width=0.35,
unit=False):
    """
    Generate a random vector for lattice constant generation. The ratios between
    x, y, and z of the returned vector correspond to the ratios between a, b,
    and c. Results in a Gaussian distribution of the natural log of the ratios.

    Args:
	minvec: the bottom-left-back minimum point which can be chosen
	maxvec: the top-right-front maximum point which can be chosen
	width: the width of the normal distribution to use when choosing values.
	    Passed to np.random.normal
	unit: whether or not to normalize the vector to determinant 1

    Returns:
	a 1x3 numpy array of floats
    """
    vec = np.array(
	[
	    np.exp(np.random.normal(scale=width)),
	    np.exp(np.random.normal(scale=width)),
	    np.exp(np.random.normal(scale=width)),
	]
    )
    if unit:
	return vec / np.linalg.norm(vec)
    else:
	return vec
				*/
