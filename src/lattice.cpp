#include "lattice.hpp"

#include <eigen3/Eigen/src/Core/Matrix.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
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

double gaussian(const double min, const double max, const double sigma) {
  double center = (max + min) * 0.5;
  double delta = fabs(max - min) * 0.5;
  double ratio = delta / sigma;
  std::normal_distribution<double> rnd{center, ratio};
  while (true) {
    double x = rnd(rand_generator);
    if (x > min && x < max) return x;
  }
  return center;
}

matrix3d parametric_to_matrix(const std::array<double, 6>& cell_params,
			      bool radians) {
  double a = cell_params[0];
  double b = cell_params[1];
  double c = cell_params[2];

  double alpha = cell_params[3];
  double beta = cell_params[4];
  double gamma = cell_params[5];
  if (!radians) {
    alpha *= xtal_consts::RAD;
    beta *= xtal_consts::RAD;
    gamma *= xtal_consts::RAD;
  }

  matrix3d matrix = matrix3d::Zero();
  double c1 = c * cos(beta);
  double c2 = (c * (cos(alpha) - (cos(beta) * cos(gamma)))) / sin(gamma);
  matrix(0, 0) = a;
  matrix(1, 0) = b * cos(gamma);
  matrix(1, 0) = b * sin(gamma);
  matrix(2, 0) = c1;
  matrix(2, 1) = c2;
  matrix(2, 2) = sqrt(c * c - c1 * c1 - c2 * c2);

  return matrix;
}

std::array<double, 6> matrix_to_parametric(const matrix3d& matrix,
					   bool radians) {
  std::array<double, 6> cell_parameters{0, 0, 0, 0, 0, 0};

  cell_parameters[0] = matrix(0, Eigen::all).norm();
  cell_parameters[1] = matrix(1, Eigen::all).norm();
  cell_parameters[2] = matrix(2, Eigen::all).norm();

  cell_parameters[3] = angle(matrix(1, Eigen::all), matrix(2, Eigen::all));
  cell_parameters[4] = angle(matrix(0, Eigen::all), matrix(2, Eigen::all));
  cell_parameters[5] = angle(matrix(0, Eigen::all), matrix(1, Eigen::all));

  if (not radians) {
    cell_parameters[3] *= xtal_consts::DEG;
    cell_parameters[4] *= xtal_consts::DEG;
    cell_parameters[5] *= xtal_consts::DEG;
  }

  return cell_parameters;
}

std::array<double, 6> generate_cellpara_0D(LatticeType ltype, double volume,
					   double maxattempts = 100) {
  double radius = cbrt((3 * volume) / (xtal_consts::PI * 4));
  double angle = 0.5 * xtal_consts::PI;
  double alpha = xtal_consts::PI / 2;
  double beta = xtal_consts::PI / 2;
  double gamma = xtal_consts::PI / 2;
  double x = (4.0 / 3.0) * xtal_consts::PI;

  switch (ltype) {
    case LatticeType::ellipsoidal:
      for (size_t i = 0; i < maxattempts; ++i) {
	auto vec = random_vector();
	double c = vec[2] / (vec[0] * vec[1]) * cbrt(volume / x);
	double a = sqrt((volume / x) / c);
	double b = sqrt((volume / x) / c);
	if ((a / c < 10.0) and (c / a < 10.0))
	  return {a, b, c, alpha, beta, gamma};
      }
    case LatticeType::spherical:
      return {radius, radius, radius, angle, angle, angle};
  }
  return {radius, radius, radius, angle, angle, angle};
}

std::array<double, 6> generate_cellpara_1D(LatticeType ltype, double volume,
					   double min_angle = xtal_consts::PI /
							      6,
					   double minvec = 1.2, double area = 0,
					   double maxattempts = 100) {
  int PA = 3;
  double max_angle = xtal_consts::PI - min_angle;
  for (size_t i = 0; i < maxattempts; ++i) {
    std::array<double, 3> abc{1, 1, 1};
    double thickness = 0;
    if (area == 0) {
      vector3d v = random_vector();
      thickness = cbrt(volume) * (v[0] / (v[0] * v[1] * v[2]));
    } else {
      thickness = volume / area;
    }
    abc[PA - 1] = thickness;
    double alpha = xtal_consts::PI / 2;
    double beta = xtal_consts::PI / 2;
    double gamma = xtal_consts::PI / 2;
    switch (ltype) {
      case LatticeType::triclinic:
	matrix3d mat = random_shear_matrix(0.2);
	auto [a, b, c, alpha, beta, gamma] = matrix_to_parametric(mat);
	double x = sqrt(1 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta) -
			cos(gamma) * cos(gamma));
	abc[PA - 1] =
	    abc[PA - 1] / x;  // scale thickness by outer product of vectors
	double ab = volume / (abc[PA - 1] * x);
	double ratio = a / b;
	switch (PA) {
	  case 3:
	    abc[0] = sqrt(ab * ratio);
	    abc[1] = sqrt(ab / ratio);
	    break;
	  case 2:
	    abc[0] = sqrt(ab * ratio);
	    abc[2] = sqrt(ab / ratio);
	    break;
	  case 1:
	    abc[1] = sqrt(ab * ratio);
	    abc[2] = sqrt(ab / ratio);
	    break;
	}
	break;
      case LatticeType::monoclinic:
	alpha = gaussian(min_angle, max_angle);
	x = sin(alpha);
	ab = volume / abc[PA - 1] * x;
	ratio = a / b;
	switch (PA) {
	  case 3:
	    abc[0] = sqrt(ab * ratio);
	    abc[1] = sqrt(ab / ratio);
	    break;
	  case 2:
	    abc[0] = sqrt(ab * ratio);
	    abc[2] = sqrt(ab / ratio);
	    break;
	  case 1:
	    abc[1] = sqrt(ab * ratio);
	    abc[2] = sqrt(ab / ratio);
	    break;
	}
	break;
      case LatticeType::orthorhombic:
	vector3d vec = random_vector();
	switch (PA) {
	  case 3:
	    ratio = abs(vec[0] / vec[1]);
	    abc[1] = sqrt(volume / (thickness * ratio));
	    abc[0] = abc[1] * ratio;
	    break;
	  case 2:
	    ratio = abs(vec[0] / vec[2]);
	    abc[2] = sqrt(volume / (thickness * ratio));
	    abc[0] = abc[2] * ratio;
	    break;
	  case 1:
	    ratio = abs(vec[1] / vec[2]);
	    abc[2] = sqrt(volume / (thickness * ratio));
	    abc[1] = abc[2] * ratio;
	    break;
	}
	break;
      case LatticeType::tetragonal:
	switch (PA) {
	  case 3:
	    abc[0] = sqrt(volume / thickness);
	    abc[1] = sqrt(volume / thickness);
	    break;
	  case 2:
	    abc[0] = abc[1];
	    abc[2] = volume / (abc[PA - 1] * abc[PA - 1]);
	    break;
	  case 1:
	    abc[1] = abc[0];
	    abc[2] = volume / (abc[PA - 1] * abc[PA - 1]);
	    break;
	}
	break;
      case LatticeType::hexagonal:
      case LatticeType::trigonal:
	gamma = xtal_consts::PI / 3 * 2;
	x = sqrt(3.0) / 2.0;
	switch (PA) {
	  case 3:
	    abc[0] = sqrt((volume / x) / abc[PA - 1]);
	    abc[1] = sqrt((volume / x) / abc[PA - 1]);
	    break;
	  case 2:
	    abc[0] = abc[1];
	    abc[2] = (volume / x) / (thickness * thickness);
	    break;
	  case 1:
	    abc[1] = abc[0];
	    abc[2] = (volume / x) / (thickness * thickness);
	    break;
	}
    }
    auto para =
	std::array<double, 6>{abc[0], abc[1], abc[2], alpha, beta, gamma};
    double a = abc[0];
    double b = abc[1];
    double c = abc[2];
    double maxvec = (a * b * c) / (minvec * minvec);
    double min_l = minvec;
    double mid_l = min_l;
    double max_l = mid_l;
    double l_min = std::min(std::min(a, b), c);
    double l_max = std::max(std::max(a, b), c);
  }
