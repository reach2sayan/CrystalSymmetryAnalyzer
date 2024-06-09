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

enum class LatticeUniqueAxis : int { a = 0, b = 1, c = 2 };
typedef LatticeUniqueAxis LUniqAx;

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

/*
enum class PeriodicBoundaryCondition : int {
  xyz = 111,
  xy = 110,
  yz = 11,
  xz = 110,
  x = 1,
  y = 2,
  z = 3
};
*/

struct LatticeParams {
  double area = 0;
  bool random = true;
  bool allow_volume_reset = false;
  LatticeUniqueAxis unique_axis = LatticeUniqueAxis::a;
  vector3d max_l{};
  vector3d min_l{};
  vector3d mid_l{};
};

class Lattice {
 public:
  Lattice(const LType _ltype, const double _volume = xtal_consts::FEMPTY,
	  const matrix3d& _matrix = {},
	  const PBC& _pbc = (PBC() << 1, 1, 1).finished(), bool _random = true,
	  bool _allow_volume_reset = false, const LUniqAx _uax = LUniqAx::c);

  Lattice(const Lattice& other) = default;
  Lattice(Lattice&& other) = default;

  int get_dimension() const {
    return std::accumulate(std::begin(pbc), std::end(pbc), 0);
  }
  int get_dofs() const;
  std::pair<vector<vector3d>, std::vector<double>> get_lengths() const;

 private:
  PBC pbc = (PBC() << 1, 1, 1).finished();
  LType ltype = LType::cubic;
  bool random = true;
  bool allow_volume_reset = false;
  int dim;
  LUniqAx unique_axis = LUniqAx::c;
  matrix3d norm_matrix;
  matrix3d stress_normalization_matrix;
  matrix3d matrix;
  std::vector<std::pair<int, int>> stress_indices;

  double a_tol;
  double volume = 0;
  int dof;
};

#endif
