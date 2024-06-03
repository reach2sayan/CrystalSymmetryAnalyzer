#ifndef __EIGEN__
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#endif

#include <unordered_map>

using vector6d = Eigen::Matrix<double, 6, 1>;
using matrix6d = Eigen::Matrix<double, 6, 6>;
using matrix3d = Eigen::Matrix<double, 3, 3>;

std::vector<std::pair<int, int>> voigt_notation = {{0, 0}, {1, 1}, {2, 2},
						   {1, 2}, {0, 2}, {0, 1}};

inline int full_3x3_to_Voight_6_index(int i, int j) {
  return i == j ? i : 6 - i - j;
}

matrix3d Voigt_6_to_full_3x3_strain(const vector6d& strain_vector) {
  double e1 = strain_vector.transpose()[0];
  double e2 = strain_vector.transpose()[1];
  double e3 = strain_vector.transpose()[2];
  double e4 = strain_vector.transpose()[3];
  double e5 = strain_vector.transpose()[4];
  double e6 = strain_vector.transpose()[5];

  matrix3d retmatrix;
  retmatrix << 1.0 + e1, 0.5 * e6, 0.5 * e5, 0.5 * e6, 1.0 + e2, 0.5 * e4,
      0.5 * e5, 0.5 * e4, 1.0 + e3;
  return retmatrix;
}

matrix3d Voigt_6_to_full_3x3_stress(const vector6d& stress_vector) {
  return stress_vector.transpose().reshaped(3, 3).eval();
}

enum class CijSymmetryTypes : int {
  CUBIC,
  TRIGONAL_H,
  TRIGONAL_L,
  TETRAGONAL_H,
  TETRAGONAL_L,
  ORTHORHOMBIC,
  MONOCLINIC,
  TRICLINIC,
  HEXAGONAL,
  NONE
};

matrix6d Cij_symmetry_cubic =
    (matrix6d() << 1, 7, 7, 0, 0, 0, 7, 1, 7, 0, 0, 0, 7, 7, 1, 0, 0, 0, 0, 0,
     0, 4, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 4)
	.finished();

matrix6d Cij_symmetry_trigonal_high =
    (matrix6d() << 1, 7, 8, 9, 10, 0, 7, 1, 8, 0, -9, 0, 8, 8, 3, 0, 0, 0, 9,
     -9, 0, 4, 0, 0, 10, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 6)
	.finished();

matrix6d Cij_symmetry_trigonal_low =
    (matrix6d() << 1, 7, 8, 9, 10, 0, 7, 1, 8, -9, -10, 0, 8, 8, 3, 0, 0, 0, 9,
     -9, 0, 4, 0, -10, 10, -10, 0, 0, 4, 9, 0, 0, 0, -10, 9, 6)
	.finished();

matrix6d Cij_symmetry_triclinic =
    (matrix6d() << 1, 7, 8, 9, 10, 11, 7, 2, 12, 13, 14, 15, 8, 12, 3, 16, 17,
     18, 9, 13, 16, 4, 19, 20, 10, 14, 17, 19, 5, 21, 11, 15, 18, 20, 21, 6)
	.finished();

matrix6d Cij_symmetry_tetragonal_high =
    (matrix6d() << 1, 7, 8, 0, 0, 0, 7, 1, 8, 0, 0, 0, 8, 8, 3, 0, 0, 0, 0, 0,
     0, 4, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 6)
	.finished();

matrix6d Cij_symmetry_tetragonal_low =
    (matrix6d() << 1, 7, 8, 0, 0, 11, 7, 1, 8, 0, 0, -11, 8, 8, 3, 0, 0, 0, 0,
     0, 0, 4, 0, 0, 0, 0, 0, 0, 4, 0, 11, -11, 0, 0, 0, 6)
	.finished();

matrix6d Cij_symmetry_orthorhombic =
    (matrix6d() << 1, 7, 8, 0, 0, 0, 7, 2, 12, 0, 0, 0, 8, 12, 3, 0, 0, 0, 0, 0,
     0, 4, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 6)
	.finished();

matrix6d Cij_symmetry_monoclinic =
    (matrix6d() << 1, 7, 8, 0, 10, 0, 7, 2, 12, 0, 14, 0, 8, 12, 3, 0, 17, 0, 0,
     0, 0, 4, 0, 20, 10, 14, 17, 0, 5, 0, 0, 0, 0, 20, 0, 6)
	.finished();

const std::unordered_map<CijSymmetryTypes, matrix6d> Cij_symmetry = {
    {CijSymmetryTypes::CUBIC, std::move(Cij_symmetry_cubic)},
    {CijSymmetryTypes::TRIGONAL_L, std::move(Cij_symmetry_trigonal_low)},
    {CijSymmetryTypes::TRIGONAL_H, std::move(Cij_symmetry_trigonal_high)},
    {CijSymmetryTypes::TETRAGONAL_H, std::move(Cij_symmetry_tetragonal_high)},
    {CijSymmetryTypes::TETRAGONAL_L, std::move(Cij_symmetry_tetragonal_low)},
    {CijSymmetryTypes::ORTHORHOMBIC, std::move(Cij_symmetry_orthorhombic)},
    {CijSymmetryTypes::MONOCLINIC, std::move(Cij_symmetry_monoclinic)},
    {CijSymmetryTypes::HEXAGONAL, std::move(Cij_symmetry_trigonal_high)},
    {CijSymmetryTypes::MONOCLINIC, std::move(Cij_symmetry_monoclinic)},
    {CijSymmetryTypes::TRICLINIC, std::move(Cij_symmetry_triclinic)}};

