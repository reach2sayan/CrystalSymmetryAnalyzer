#ifndef __SYMMETRY_ANALYZER_CONSTANTS_HPP__
#define __SYMMETRY_ANALYZER_CONSTANTS_HPP__

#include <array>
#include <cmath>
#include <set>

#ifndef __EIGEN__
#include <eigen3/Eigen/Dense>
#endif

namespace CrystalSymmetryConstants {
constexpr double DEFAULT_SYMMETRY_TOLERANCE = 0.01;
constexpr double DEFAULT_EQUIVALENCE_TOLERANCE = 0.001;
constexpr double PI = 3.14159265359;
constexpr double RAD = PI / 180;
constexpr double SQRT_INV_PI = 0.56418958354;
constexpr double SQRT_3_4PI = 0.4886025119;
constexpr double SQRT_15_PI = 2.18509686118;
constexpr double SQRT_5_PI = 1.26156626101;
constexpr double DEG = 180 / PI;
constexpr short VERBOSITY = 1;

Eigen::Matrix3d HEXAGONAL_CELL((Eigen::Matrix3d() << 1, -0.5, 0, 0,
				std::sqrt(3) / 2, 0, 0, 0, 1)
				   .finished());

const std::set<std::array<int, 3>> all_sym_directions = {
    std::array<int, 3>{1, 0, 0},   std::array<int, 3>{0, 1, 0},
    std::array<int, 3>{0, 0, 1},   std::array<int, 3>{1, 1, 1},
    std::array<int, 3>{1, -1, -1}, std::array<int, 3>{-1, 1, -1},
    std::array<int, 3>{1, -1, 0},  std::array<int, 3>{1, 1, 0},
    std::array<int, 3>{0, 1, -1},  std::array<int, 3>{0, 1, 1},
    std::array<int, 3>{-1, 0, 1},  std::array<int, 3>{1, 0, 1},
    std::array<int, 3>{1, -2, 0},  std::array<int, 3>{2, -1, 0}};

enum class coordinate_labels : int { x = 0, y = 1, z = 3 };
}  // namespace CrystalSymmetryConstants
namespace xtal_consts = CrystalSymmetryConstants;
#endif
