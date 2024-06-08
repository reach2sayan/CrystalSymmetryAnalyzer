#ifndef __SYMMETRY_ANALYZER_SPHERICAL_HARMONICS_HPP__
#define __SYMMETRY_ANALYZER_SPHERICAL_HARMONICS_HPP__

template <unsigned int N>
struct Factorial {
  static const unsigned long long value = N * Factorial<N - 1>::value;
};

// Base case specialization
template <>
struct Factorial<0> {
  static const unsigned long long value = 1;
};

constexpr const unsigned long long FACTORIAL_0 = 1;
constexpr const unsigned long long FACTORIAL_1 = 1;
constexpr const unsigned long long FACTORIAL_2 = 2;
constexpr const unsigned long long FACTORIAL_3 = Factorial<3>::value;
constexpr const unsigned long long FACTORIAL_4 = Factorial<4>::value;
constexpr const unsigned long long FACTORIAL_5 = Factorial<5>::value;
constexpr const unsigned long long FACTORIAL_6 = Factorial<6>::value;
constexpr const unsigned long long FACTORIAL_7 = Factorial<7>::value;
constexpr const unsigned long long FACTORIAL_8 = Factorial<8>::value;
constexpr const unsigned long long FACTORIAL_9 = Factorial<9>::value;
constexpr const unsigned long long FACTORIAL_10 = Factorial<9>::value;
constexpr const unsigned long long FACTORIAL_11 = Factorial<11>::value;
constexpr const unsigned long long FACTORIAL_12 = Factorial<12>::value;
constexpr const unsigned long long FACTORIAL_13 = Factorial<13>::value;
constexpr const unsigned long long FACTORIAL_14 = Factorial<14>::value;
constexpr const unsigned long long FACTORIAL_15 = Factorial<15>::value;
constexpr const unsigned long long FACTORIAL_16 = Factorial<16>::value;

constexpr int FACTORIAL_CACHE_SIZE = 17;
constexpr const unsigned long long FACTORIAL_CACHE[17] = {
    FACTORIAL_0,  FACTORIAL_1,	FACTORIAL_2,  FACTORIAL_3,  FACTORIAL_4,
    FACTORIAL_5,  FACTORIAL_6,	FACTORIAL_7,  FACTORIAL_8,  FACTORIAL_9,
    FACTORIAL_10, FACTORIAL_11, FACTORIAL_12, FACTORIAL_13, FACTORIAL_14,
    FACTORIAL_15, FACTORIAL_16};

constexpr unsigned long long factorial(int val);

double spherical_harmonics(int l, int m, double phi, double theta);

#endif
