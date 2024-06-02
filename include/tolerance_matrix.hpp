#ifndef __TOLERANCE_MATRIX_HPP__
#define __TOLERANCE_MATRIX_HPP__

enum class Prototype : int { atomic = 0, molecular = 1, vdW = 2, metallic = 3 };

class ToleranceMatrix {
 private:
  double f;
  Prototype prototype;

 public:
	ToleranceMatrix(const std::vector<
};

#endif
