#ifndef __LATTICE_HPP__
#define __LATTICE_HPP__

#include <array>
#include <numeric>
#include <vector>
#define eigen_assert(X)                     \
  do {                                      \
    if (!(X)) throw std::runtime_error(#X); \
  } while (false);
#if defined(WIN32) || defined(_WIN32) || \
    defined(__WIN32) && !defined(__CYGWIN__)

#else
#ifndef __EIGEN__
#include <eigen3/Eigen/Dense>
template <typename T>
using vector = std::vector<T, Eigen::aligned_allocator<T>>;
using PBC = Eigen::Array<bool, 3, 1>;
using vector3d = Eigen::Matrix<double, 3, 1>;
using matrix3d = Eigen::Matrix<double, 3, 3>;
#endif
#endif

#ifndef __SYMMETRY_ANALYZER_CONSTANTS_HPP__
#include "constants.hpp"
#endif

#include <random>
std::random_device rand_device{};
std::mt19937 rand_generator{rand_device()};

enum class LatticeUniqueAxis : int { a = 0, b = 1, c = 2 };
typedef LatticeUniqueAxis LUniqAx;

using CellParams = std::tuple<double, double, double, double, double, double>;
const std::unordered_map<std::string, unsigned short> cell_param_index_map{
    {"a", 0}, {"b", 1}, {"c", 2}, {"alpha", 3}, {"beta", 4}, {"gamma", 5}};

enum class LatticeType : int {
  triclinic = 0,
  monoclinic = 1,
  orthorhombic = 2,
  tetragonal = 3,
  trigonal = 4,
  hexagonal = 5,
  cubic = 6,
  spherical = 100,
  ellipsoidal = 200
};
typedef LatticeType LType;

class Lattice {
 public:
  Lattice(const LType _ltype, const double _volume = xtal_consts::FEMPTY,
	  const matrix3d& _matrix = {}, const PBC& _pbc = {true, true, true},
	  bool _random = true, bool _allow_volume_reset = false,
	  const LUniqAx _uax = LUniqAx::c);

  Lattice() = delete;
  Lattice(const Lattice& other) = default;
  Lattice(Lattice&& other) = default;

  int get_dimension() const {
    return std::accumulate(std::begin(pbc), std::end(pbc), 0);
  }
  int get_dofs() const;
  std::pair<vector<vector3d>, std::vector<double>> get_lengths() const;
  vector<matrix3d> get_permutation_matrices() const;
  vector<matrix3d> get_transformation_matrices() const;
  double get_worst_angle() const;

  const matrix3d& get_matrix() const { return matrix; }
  CellParams get_parameters(bool degree = false) const;

 private:
  PBC pbc{true, true, true};
  LType ltype = LType::cubic;
  double a = 1.0;
  double b = 1.0;
  double c = 1.0;
  double alpha = 90.;
  double gamma = 90.;
  double beta = 90.;
  bool random = true;
  bool allow_volume_reset = false;
  int dim = 3;
  LUniqAx unique_axis = LUniqAx::c;
  matrix3d norm_matrix;
  matrix3d stress_normalization_matrix;
  matrix3d matrix;
  std::vector<std::pair<int, int>> stress_indices;

  double a_tol = xtal_consts::DEFAULT_EQUIVALENCE_TOLERANCE;
  double volume = 1.0;
  int dof = 3;
};

CellParams generate_cell_parameters(const LatticeType ltype, const double vol,
				    const double minvec = 1.2,
				    const double minangle = xtal_consts::PI / 6,
				    const double max_ratio = 10.0,
				    const int max_attempts = 100);

matrix3d random_shear_matrix(double width = 1.0, bool unitary = false);
vector3d random_vector(const vector3d& minvector = {0.0, 0.0, 0.0},
		       const vector3d& maxvector = {1.0, 1.0, 1.0},
		       double width = 0.35, bool unit = false);

double gaussian(const double min, const double max, const double sigma = 3.0);

matrix3d parametric_to_matrix(const std::array<double, 6>& cell_params,
			      bool radians = true);

std::array<double, 6> matrix_to_parametric(const matrix3d& matrix,
					   bool radians = true);

#endif
