#ifndef __SYMMETRY_ANALYZER_SPHERICAL_HARMONICS_HPP__
#define __SYMMETRY_ANALYZER_SPHERICAL_HARMONICS_HPP__

template <long long N>
struct Factorial {
  enum { value = N * Factorial<N - 1>::value };
};

template <>
struct Factorial<0> {
  enum { value = 1 };
};

constexpr double FACTORIAL_0 = 1;
constexpr double FACTORIAL_1 = 1;
constexpr double FACTORIAL_2 = 2;
constexpr double FACTORIAL_3 = Factorial<3>::value;
constexpr double FACTORIAL_4 = Factorial<4>::value;
constexpr double FACTORIAL_5 = Factorial<5>::value;
constexpr double FACTORIAL_6 = Factorial<6>::value;
constexpr double FACTORIAL_7 = Factorial<7>::value;
constexpr double FACTORIAL_8 = Factorial<8>::value;
constexpr double FACTORIAL_9 = Factorial<9>::value;
constexpr double FACTORIAL_10 = Factorial<9>::value;
constexpr double FACTORIAL_11 = Factorial<11>::value;
constexpr double FACTORIAL_12 = Factorial<12>::value;
constexpr double FACTORIAL_13 = Factorial<13>::value;
constexpr double FACTORIAL_14 = Factorial<14>::value;
constexpr double FACTORIAL_15 = Factorial<15>::value;
constexpr double FACTORIAL_16 = Factorial<16>::value;

constexpr int FACTORIAL_CACHE_SIZE = 17;
constexpr const double FACTORIAL_CACHE[17] = {
    FACTORIAL_0,  FACTORIAL_1,	FACTORIAL_2,  FACTORIAL_3,  FACTORIAL_4,
    FACTORIAL_5,  FACTORIAL_6,	FACTORIAL_7,  FACTORIAL_8,  FACTORIAL_9,
    FACTORIAL_10, FACTORIAL_11, FACTORIAL_12, FACTORIAL_13, FACTORIAL_14,
    FACTORIAL_15, FACTORIAL_16};
constexpr double factorial(int val);

double spherical_harmonics(int l, int m, double phi, double theta);

#endif
