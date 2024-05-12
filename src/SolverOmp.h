#ifndef SOLVEROMP_H_
#define SOLVEROMP_H_

#include <cstddef>
#include <fstream>
#include <vector>
class Solver {
 public:
  Solver(const std::ptrdiff_t nx, const std::ptrdiff_t nt, const double c,
         const std::ptrdiff_t interval);

  ~Solver();

  void Initialize(double (*pInitU)(const double),
                  double (*pInitUT)(const double, const double));

  void Solve();

 private:
  int mNext = 0;
  int mCurr = 1;
  int mPrev = 2;

  const std::ptrdiff_t mNx;
  const std::ptrdiff_t mNt;
  const double mC;
  const std::ptrdiff_t mInterval;
  const std::ptrdiff_t mR;

  double mDx;
  double mDt;
  double mS;
  double mS2;

  std::vector<std::vector<double>> mU;

  std::ofstream mOutputFile;
  std::ofstream mExactFile;
  std::ofstream mDiffFile;

  void Step();
  void SolveOneStep();
  /** Computes the exact solution of u(x,t) */
  std::vector<double> GetExactSolution(const double time);
  void WriteOutput(const double time, const std::vector<double>& rU,
                   const std::vector<double>& rExact,
                   const std::vector<double>& rDiff);
  /** Computes the difference between two arrays */
  static std::vector<double> GetDifference(const std::vector<double>& rArr1,
                                           const std::vector<double>& rArr2);
  /**
   * Writes array to file.
   *
   * @param time  Current time.
   * @param rArr  Array containing values of u.
   * @parem rOutfile  Output file.
   */
  static void WriteToFile(const double time, const std::vector<double>& rArr,
                          std::ofstream& rOutfile);
};

#endif  // SOLVEROMP_H_
