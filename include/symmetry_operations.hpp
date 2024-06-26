#ifndef __SYMMETRY_ANALYZER_OPERATIONS_HPP__
#define __SYMMETRY_ANALYZER_OPERATIONS_HPP__

#include <variant>
#include <vector>

#ifndef __EIGEN__
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
template <typename T>
using vector = std::vector<T, Eigen::aligned_allocator<T>>;
using PBC = Eigen::Array<bool, 3, 1>;
using matrix4d = Eigen::Matrix<double, 4, 4>;
using vector4d = Eigen::Matrix<double, 4, 1>;
using matrix3d = Eigen::Matrix<double, 3, 3>;
using vector3d = Eigen::Matrix<double, 3, 1>;
#endif

#ifndef __SYMMETRY_ANALYZER_CONSTANTS_HPP__
#include "constants.hpp"
#endif

class SymmetryOperations {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  SymmetryOperations(const matrix4d& _matrix, double _tol)
      : affine_transformation_matrix{_matrix}, tolerance{_tol} {}

  bool operator==(const SymmetryOperations& other);
  SymmetryOperations operator*(const SymmetryOperations& rhs);

  vector3d operator()(const vector3d& point) const;
  vector<vector3d> operator()(const vector<vector3d>& points) const;

  vector3d apply_rotation_only(const vector3d& v) const;

  inline matrix3d rotation_matrix() const {
    return affine_transformation_matrix.block(0, 0, 2, 2);
  }
  inline vector3d translation_vector() const {
    return affine_transformation_matrix.block(0, 3, 2, 3);
  }

  SymmetryOperations inverse() const {
    return {affine_transformation_matrix.inverse(), tolerance};
  }

  bool are_symmetrically_related(
      const vector3d& point_a, const vector3d& point_b,
      double tol = xtal_consts::DEFAULT_EQUIVALENCE_TOLERANCE) const;

  template <size_t rank>
  Eigen::Tensor<double, rank> transform_tensor(
      const Eigen::Tensor<double, rank>& tensor) const;

  std::pair<bool, bool> are_symmetrically_related_vectors(
      const vector3d& from_a, const vector3d& to_a, const vector3d& r_a,
      const vector3d& from_b, const vector3d& to_b, const vector3d& r_b,
      double tol = xtal_consts::DEFAULT_EQUIVALENCE_TOLERANCE);

  const matrix4d& get_matrix() const { return affine_transformation_matrix; }

 private:
  matrix4d affine_transformation_matrix;
  double tolerance = xtal_consts::DEFAULT_SYMMETRY_TOLERANCE;
};

SymmetryOperations from_rotation_and_translation(
    const matrix3d& rotation_matrix, const vector3d& translation_matrix,
    double tolerance = xtal_consts::DEFAULT_SYMMETRY_TOLERANCE);

SymmetryOperations from_axis_angle_and_translation(
    const vector3d axis, double angle, bool angle_in_radians = false,
    vector3d translation_vector = vector3d::Zero());

SymmetryOperations from_origin_axis_angle(const vector3d& origin,
					  const vector3d& axis, double angle,
					  bool angle_in_radians = false);

SymmetryOperations reflection(const vector3d& normal,
			      const vector3d& origin = vector3d::Zero());

SymmetryOperations rotoreflection(const vector3d& axis, double angle,
				  const vector3d& origin = vector3d::Zero());

#endif
