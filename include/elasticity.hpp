#ifndef __SYMMETRY_ANALYZER_ELASTICITY_HPP__
#define __SYMMETRY_ANALYZER_ELASTICITY_HPP__

#ifndef __EIGEN__
#include <array>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
using vector6d = Eigen::Matrix<double, 6, 1>;
using vector3d = Eigen::Matrix<double, 3, 1>;
using matrix6d = Eigen::Matrix<double, 6, 6>;
using matrix3d = Eigen::Matrix<double, 3, 3>;
using tensor4r = Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3, 3>>;
#endif
#include <set>
#include <unordered_map>

const std::set<std::pair<int, int>> voigt_notation = {{0, 0}, {1, 1}, {2, 2},
						      {1, 2}, {0, 2}, {0, 1}};

inline int full_3x3_to_Voight_6_index(int i, int j) {
  assert(i < 3 && j < 3);
  return i == j ? i : 6 - i - j;
}

matrix3d Voigt_6_to_full_3x3_strain(const vector6d& strain_vector);
matrix3d Voigt_6_to_full_3x3_stress(const vector6d& stress_vector);
vector6d full_3x3_to_Voigt_6_strain(const matrix3d& strain_matrix);
vector6d full_3x3_to_Voigt_6_stress(const matrix3d& stress_matrix);
matrix6d full_3x3x3x3_to_Voigt_6x6(const tensor4r& C);
vector3d Voigt_6x6_to_cubic(const matrix6d& C);

using full_3x3_to_Voigt_6 = vector6d (*)(const matrix3d&);
std::tuple<double, double, double> __invariants_impl(const vector6d& voigt);
std::tuple<double, double, double> __invariants_impl(double sxx, double syy,
						     double szz, double syz,
						     double sxz, double sxy);
std::tuple<double, double, double> __invariants_impl(const matrix3d& matrix,
						     full_3x3_to_Voigt_6);

template <typename... Args>
constexpr std::tuple<double, double, double> invariants(Args... args);

matrix6d rotate_cubic_elastic_constants(double C11, double C12, double C44,
					const matrix3d& A, double tol = 1e-6);

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

const matrix6d Cij_symmetry_cubic =
    (matrix6d() << 1, 7, 7, 0, 0, 0, 7, 1, 7, 0, 0, 0, 7, 7, 1, 0, 0, 0, 0, 0,
     0, 4, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 4)
	.finished();

const matrix6d Cij_symmetry_trigonal_high =
    (matrix6d() << 1, 7, 8, 9, 10, 0, 7, 1, 8, 0, -9, 0, 8, 8, 3, 0, 0, 0, 9,
     -9, 0, 4, 0, 0, 10, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 6)
	.finished();

const matrix6d Cij_symmetry_trigonal_low =
    (matrix6d() << 1, 7, 8, 9, 10, 0, 7, 1, 8, -9, -10, 0, 8, 8, 3, 0, 0, 0, 9,
     -9, 0, 4, 0, -10, 10, -10, 0, 0, 4, 9, 0, 0, 0, -10, 9, 6)
	.finished();

const matrix6d Cij_symmetry_triclinic =
    (matrix6d() << 1, 7, 8, 9, 10, 11, 7, 2, 12, 13, 14, 15, 8, 12, 3, 16, 17,
     18, 9, 13, 16, 4, 19, 20, 10, 14, 17, 19, 5, 21, 11, 15, 18, 20, 21, 6)
	.finished();

const matrix6d Cij_symmetry_tetragonal_high =
    (matrix6d() << 1, 7, 8, 0, 0, 0, 7, 1, 8, 0, 0, 0, 8, 8, 3, 0, 0, 0, 0, 0,
     0, 4, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 6)
	.finished();

const matrix6d Cij_symmetry_tetragonal_low =
    (matrix6d() << 1, 7, 8, 0, 0, 11, 7, 1, 8, 0, 0, -11, 8, 8, 3, 0, 0, 0, 0,
     0, 0, 4, 0, 0, 0, 0, 0, 0, 4, 0, 11, -11, 0, 0, 0, 6)
	.finished();

const matrix6d Cij_symmetry_orthorhombic =
    (matrix6d() << 1, 7, 8, 0, 0, 0, 7, 2, 12, 0, 0, 0, 8, 12, 3, 0, 0, 0, 0, 0,
     0, 4, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 6)
	.finished();

const matrix6d Cij_symmetry_monoclinic =
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
#endif
