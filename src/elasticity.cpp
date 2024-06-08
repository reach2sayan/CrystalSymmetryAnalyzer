#include "elasticity.hpp"

#include <algorithm>
#include <cassert>
#include <vector>

matrix3d Voigt_6_to_full_3x3_strain(const vector6d& strain_vector) {
  double e1 = strain_vector.transpose()[0];
  double e2 = strain_vector.transpose()[1];
  double e3 = strain_vector.transpose()[2];
  double e4 = strain_vector.transpose()[3];
  double e5 = strain_vector.transpose()[4];
  double e6 = strain_vector.transpose()[5];

  return (matrix3d() << 1.0 + e1, 0.5 * e6, 0.5 * e5, 0.5 * e6, 1.0 + e2,
	  0.5 * e4, 0.5 * e5, 0.5 * e4, 1.0 + e3)
      .finished();
}

matrix3d Voigt_6_to_full_3x3_stress(const vector6d& stress_vector) {
  double s1 = stress_vector.transpose()[0];
  double s2 = stress_vector.transpose()[1];
  double s3 = stress_vector.transpose()[2];
  double s4 = stress_vector.transpose()[3];
  double s5 = stress_vector.transpose()[4];
  double s6 = stress_vector.transpose()[5];

  return (matrix3d() << s1, s6, s5, s6, s2, s4, s5, s4, s3)
      .finished()
      .transpose();
}

vector6d full_3x3_to_Voigt_6_strain(const matrix3d& strain_matrix) {
  assert(strain_matrix.transpose() == strain_matrix);
  return (vector6d() << strain_matrix(0, 0) - 1.0, strain_matrix(1, 1) - 1.0,
	  strain_matrix(2, 2) - 1.0, strain_matrix(1, 2) + strain_matrix(2, 1),
	  strain_matrix(0, 2) + strain_matrix(2, 0),
	  strain_matrix(0, 1) + strain_matrix(1, 0))
      .finished();
}

vector6d full_3x3_to_Voigt_6_stress(const matrix3d& stress_matrix) {
  assert(stress_matrix.transpose() == stress_matrix);
  return (vector6d() << stress_matrix(0, 0), stress_matrix(1, 1),
	  stress_matrix(2, 2),
	  (stress_matrix(1, 2) + stress_matrix(2, 1)) / 2.0,
	  (stress_matrix(0, 2) + stress_matrix(2, 0)) / 2.0,
	  (stress_matrix(0, 1) + stress_matrix(1, 0)) / 2.0)
      .finished();
}

matrix6d full_3x3x3x3_to_Voigt_6x6(const tensor4r& C) {
  double tol = 1e-3;
  matrix6d voigt = matrix6d::Zero();

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      auto [k, l] = voigt_notation[i];
      auto [m, n] = voigt_notation[j];
      voigt(i, j) = C(k, l, m, n);

      assert(fabs(voigt(i, j) - C(m, n, k, l)) < tol);
      assert(fabs(voigt(i, j) - C(l, k, m, n)) < tol);
      assert(fabs(voigt(i, j) - C(k, l, n, m)) < tol);
      assert(fabs(voigt(i, j) - C(m, n, l, k)) < tol);
      assert(fabs(voigt(i, j) - C(n, m, k, l)) < tol);
      assert(fabs(voigt(i, j) - C(l, k, n, m)) < tol);
      assert(fabs(voigt(i, j) - C(n, m, l, k)) < tol);
    }
  }
  return voigt;
}

vector3d Voigt_6x6_to_cubic(const matrix6d& C) {
  double tol = 1e-6;
  vector3d C11s = (vector3d() << C(0, 0), C(1, 1), C(2, 2)).finished();
  vector3d C12s = (vector3d() << C(1, 2), C(0, 2), C(0, 1)).finished();
  vector3d C44s = (vector3d() << C(3, 3), C(4, 4), C(5, 5)).finished();

  matrix6d C_check = matrix6d::Zero();
  C_check.diagonal() = C.diagonal();
  C_check(Eigen::seq(0, 3), Eigen::seq(0, 3)) =
      C(Eigen::seq(0, 3), Eigen::seq(0, 3));

  matrix6d diff = C - C_check;
  assert(std::none_of(diff.reshaped().begin(), diff.reshaped().end(),
		      [tol](auto diff) { return fabs(diff) > tol; }));

  double C11 = C11s.mean();
  double C12 = C12s.mean();
  double C44 = C44s.mean();

  assert(std::none_of(C11s.begin(), C11s.end(), [C11, tol](double val) {
    return fabs(val - C11) > tol;
  }));
  assert(std::none_of(C12s.begin(), C12s.end(), [C12, tol](double val) {
    return fabs(val - C12) > tol;
  }));
  assert(std::none_of(C44s.begin(), C44s.end(), [C44, tol](double val) {
    return fabs(val - C44) > tol;
  }));

  return (vector3d() << C11, C12, C44).finished();
}

matrix6d cubic_to_Voigt_6x6(double C11, double C12, double C44) {
  return (matrix6d() << C11, C12, C12, 0, 0, 0, C12, C11, C12, 0, 0, 0, C12,
	  C12, C11, 0, 0, 0, 0, 0, 0, C44, 0, 0, 0, 0, 0, 0, C44, 0, 0, 0, 0, 0,
	  0, C44)
      .finished();
}

std::tuple<double, double, double> __invariants_impl(const vector6d& voigt) {
  double I1 = voigt[0] + voigt[1] + voigt[2];
  double I2 = voigt[0] * voigt[1] + voigt[1] * voigt[2] + voigt[2] * voigt[0] -
	      voigt[3] * voigt[3] - voigt[4] * voigt[4] - voigt[5] * voigt[5];
  double I3 = voigt[0] * voigt[1] * voigt[2] +
	      2 * voigt[3] * voigt[4] * voigt[5] -
	      voigt[3] * voigt[3] * voigt[2] - voigt[4] * voigt[4] * voigt[0] -
	      voigt[5] * voigt[5] * voigt[1];
  return {I1, I2, I3};
}

std::tuple<double, double, double> __invariants_impl(double sxx, double syy,
						     double szz, double syz,
						     double sxz, double sxy) {
  const vector6d voigt =
      (vector6d() << sxx, syy, szz, syz, sxz, sxy).finished();

  return __invariants_impl(voigt);
}

std::tuple<double, double, double> __invariants_impl(const matrix3d& matrix,
						     full_3x3_to_Voigt_6 func) {
  vector6d voigt = func(matrix);
  return __invariants_impl(voigt);
}

template <typename... Args>
constexpr std::tuple<double, double, double> invariants(Args... args) {
  auto [I1, I2, I3] = __invariants_impl(std::forward<Args>(args)...);

  double J1 = -I1 / 3;
  double J2 = I1 * I1 / 3 - I2;
  double J3 = 2 * I1 * I1 * I1 / 27 - I1 * I2 / 3 + I3;

  return {-J1, sqrt(2 * J2 / 3), J3};
}

// TODO : Incomplete
matrix6d rotate_cubic_elastic_constants(double C11, double C12, double C44,
					matrix3d A, double tol = 1e-6) {
  vector<double> C;
  for (const auto [i, j] : voigt_notation) {
    for (const auto [k, l] : voigt_notation) {
      double h = 0;
      if (i == j && k == l) h += C12;  // la
      if (i == k && j == l) h += C44;  // mu
      if (i == l && j == k) h += C44;
      h += (C11 - C12 - 2 * C44);  //*np.sum(A[i,:]*A[j,:]*A[k,:]*A[l,:])
      C.push_back(h);
    }
  }
  return Eigen::Map<matrix6d>(C.data());
}
