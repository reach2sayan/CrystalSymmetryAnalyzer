#if !defined(__SYMMETRY_ANALYZER_CONSTANTS_HPP__)
#include "constants.hpp"
#endif

#if !defined(__SYMMETRY_ANALYZER_OPERATIONS_HPP__)
#include "symmetry_operations.hpp"
#endif

#include <algorithm>
#include <cmath>

template <typename T>
constexpr T to_radians(T num) {
  return num * xtal_consts::PI / 180;
}

bool SymmetryOperations::operator==(const SymmetryOperations& other) {
  if (affine_transformation_matrix.size() !=
      other.affine_transformation_matrix.size())
    return false;

  const auto this_reshaped = affine_transformation_matrix.reshaped();
  const auto other_reshaped = affine_transformation_matrix.reshaped();
  size_t size =
      affine_transformation_matrix.rows() + affine_transformation_matrix.cols();
  for (size_t i = 0; i < size; i++)
    if (std::fabs(this_reshaped[i] - other_reshaped[i]) > tolerance)
      return false;
  return true;
}

SymmetryOperations SymmetryOperations::operator*(
    const SymmetryOperations& rhs) {
  return {affine_transformation_matrix * rhs.affine_transformation_matrix,
	  std::min(tolerance, rhs.tolerance)};
}

SymmetryOperations from_axis_angle_and_translation(
    const vector3d axis, double angle, bool angle_in_radians,
    vector3d translation_vector) {
  double ang = angle_in_radians ? angle : to_radians(angle);
  double cos_a = std::cos(ang);
  double sin_a = std::sin(ang);
  vector3d unit_vec = axis / axis.norm();

  matrix3d rot_mat =
      (matrix3d() << cos_a + unit_vec[0] * unit_vec[0] * (1 - cos_a),
       unit_vec[0] * unit_vec[1] * (1 - cos_a) - unit_vec[2] * sin_a,
       unit_vec[0] * unit_vec[2] * (1 - cos_a) + unit_vec[1] * sin_a,
       unit_vec[0] * unit_vec[1] * (1 - cos_a) + unit_vec[2] * sin_a,
       cos_a + unit_vec[1] * unit_vec[1] * (1 - cos_a),
       unit_vec[1] * unit_vec[2] * (1 - cos_a) - unit_vec[0] * sin_a,
       unit_vec[0] * unit_vec[2] * (1 - cos_a) - unit_vec[1] * sin_a,
       unit_vec[1] * unit_vec[2] * (1 - cos_a) + unit_vec[0] * sin_a,
       cos_a + unit_vec[2] * unit_vec[2] * (1 - cos_a))
	  .finished();

  return from_rotation_and_translation(rot_mat, translation_vector);
}

SymmetryOperations from_origin_axis_angle(const vector3d& origin,
					  const vector3d& axis, double angle,
					  bool angle_in_radians) {
  double theta = angle_in_radians ? angle : to_radians(angle);

  double a = origin(0);
  double b = origin(1);
  double c = origin(2);

  double ax_u = axis(0);
  double ax_v = axis(1);
  double ax_w = axis(2);

  double u2 = ax_u * ax_u;
  double v2 = ax_v * ax_v;
  double w2 = ax_w * ax_w;

  double cos_t = std::cos(theta);
  double sin_t = std::sin(theta);

  double l2 = u2 + v2 + w2;
  double lsqrt = std::sqrt(l2);

  matrix4d affine_mat =
      (matrix4d() << (u2 + (v2 + w2) * cos_t) / l2,
       (ax_u * ax_v * (1 - cos_t) - ax_w * lsqrt * sin_t) / l2,
       (ax_u * ax_w * (1 - cos_t) + ax_v * lsqrt * sin_t) / l2,
       (a * (v2 + w2) - ax_u * (b * ax_v + c * ax_w) +
	(ax_u * (b * ax_v + c * ax_w) - a * (v2 + w2)) * cos_t +
	(b * ax_w - c * ax_v) * lsqrt * sin_t) /
	   l2,

       (ax_u * ax_v * (1 - cos_t) + ax_w * lsqrt * sin_t) / l2,
       (v2 + (u2 + w2) * cos_t) / l2,
       (ax_v * ax_w * (1 - cos_t) - ax_u * lsqrt * sin_t) / l2,
       (b * (u2 + w2) - ax_v * (a * ax_u + c * ax_w) +
	(ax_v * (a * ax_u + c * ax_w) - b * (u2 + w2)) * cos_t +
	(c * ax_u - a * ax_w) * lsqrt * sin_t) /
	   l2,

       (ax_u * ax_w * (1 - cos_t) - ax_v * lsqrt * sin_t) / l2,
       (ax_v * ax_w * (1 - cos_t) + ax_u * lsqrt * sin_t) / l2,
       (w2 + (u2 + v2) * cos_t) / l2,
       (c * (u2 + v2) - ax_w * (a * ax_u + b * ax_v) +
	(ax_w * (a * ax_u + b * ax_v) - c * (u2 + v2)) * cos_t +
	(a * ax_v - b * ax_u) * lsqrt * sin_t) /
	   l2,

       0, 0, 0, 1)
	  .finished();

  return {affine_mat, xtal_consts::DEFAULT_SYMMETRY_TOLERANCE};
}

SymmetryOperations from_rotation_and_translation(
    const matrix3d& rotation_matrix, const vector3d& translation_matrix,
    double tolerance) {
  matrix4d affine_matrix = matrix4d::Identity();
  affine_matrix.block(0, 0, 2, 2) = rotation_matrix;
  affine_matrix(Eigen::seq(0, 3), Eigen::seq(2, 3)) = translation_matrix;

  return {affine_matrix, tolerance};
}

vector3d SymmetryOperations::operator()(const vector3d& point) const {
  vector4d affine_point = vector4d::Ones();
  affine_point.head(3) = point;
  return (affine_transformation_matrix * affine_point).head(3);
}

vector<vector3d> SymmetryOperations::operator()(
    const vector<vector3d>& points) const {
  vector<vector3d> retvector;
  retvector.reserve(points.size());

  for (const vector3d& point : points)
    retvector.emplace_back(operator()(point));

  return retvector;
}

vector3d SymmetryOperations::apply_rotation_only(const vector3d& v) const {
  return affine_transformation_matrix.block(0, 0, 3, 3) * v;
}

bool SymmetryOperations::are_symmetrically_related(const vector3d& point_a,
						   const vector3d& point_b,
						   double tol) const {
  vector3d diff_ab = operator()(point_a) - point_b;
  vector3d diff_ba = operator()(point_b) - point_a;
  auto equal_within_tolerance = [tol](const double a) -> bool {
    return std::fabs(a) < tol;
  };
  return std::all_of(diff_ab.begin(), diff_ab.end(), equal_within_tolerance) &&
	 std::all_of(diff_ba.begin(), diff_ba.end(), equal_within_tolerance);
}

std::pair<bool, bool> SymmetryOperations::are_symmetrically_related_vectors(
    const vector3d& from_a, const vector3d& to_a, const vector3d& r_a,
    const vector3d& from_b, const vector3d& to_b, const vector3d& r_b,
    double tol) {
  vector3d from_c = operator()(from_a);
  vector3d to_c = operator()(to_a);
  return {false, false};
}

SymmetryOperations reflection(const vector3d& normal, const vector3d& origin) {
  vector3d normalized = normal.normalized();
  double u = normalized[0];
  double v = normalized[1];
  double w = normalized[2];

  matrix4d translation = matrix4d::Identity();
  translation(Eigen::seq(0, 3), Eigen::seq(0, 2)) = -1 * origin;
  double xx = 1 - 2 * u * u;
  double yy = 1 - 2 * v * v;
  double zz = 1 - 2 * w * w;
  double xy = -2 * u * v;
  double xz = -2 * u * w;
  double yz = -2 * v * w;

  matrix4d mirror_mat =
      (matrix4d() << xx, xy, xz, 0, xy, yy, yz, 0, xz, yz, zz, 0, 0, 0, 0, 1)
	  .finished();

  if (origin.norm() > 1e-06)
    mirror_mat = translation.inverse() * mirror_mat * translation;

  return {mirror_mat, xtal_consts::DEFAULT_SYMMETRY_TOLERANCE};
}

SymmetryOperations rotoreflection(const vector3d& axis, double angle,
				  const vector3d& origin) {
  SymmetryOperations rot = from_origin_axis_angle(origin, axis, angle);
  SymmetryOperations refl = reflection(axis, origin);

  return {rot.get_matrix() * refl.get_matrix(),
	  xtal_consts::DEFAULT_SYMMETRY_TOLERANCE};
}

template <size_t rank>
Eigen::Tensor<double, rank> SymmetryOperations::transform_tensor(
    const Eigen::Tensor<double, rank>& tensor) const {
  return affine_transformation_matrix.transpose() * tensor *
	 affine_transformation_matrix;
}

