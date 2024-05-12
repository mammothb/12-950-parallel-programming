#include "SolverOmp.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <utility>

#include "omp.h"

Solver::Solver(const std::ptrdiff_t nx, const std::ptrdiff_t nt, const double c,
               const std::ptrdiff_t interval)
    : mNx(nx), mNt(nt), mC(c), mInterval(interval), mR(mNx - 1) {
  mDx = 1.0 / static_cast<double>(mNx - 1);
  mDt = 1.0 / static_cast<double>(mNt - 1);
  mS = mC * mDt / mDx;
  mS2 = mS * mS;

  std::cout << "CFL condition: " << mS << std::endl;
  assert(mS < 1.0);

  for (int i = 0; i < 3; ++i) {
    mU.emplace_back(mNx);
  }

  mOutputFile.open("result.csv", std::ios::out);
  mOutputFile << std::fixed << std::setprecision(6);
  mExactFile.open("exact.csv", std::ios::out);
  mExactFile << std::fixed << std::setprecision(6);
  mDiffFile.open("diff.csv", std::ios::out);
  mDiffFile << std::fixed << std::setprecision(6);
}

Solver::~Solver() {
  mOutputFile.close();
  mExactFile.close();
  mDiffFile.close();
}

void Solver::Initialize(double (*pInitU)(const double),
                        double (*pInitUT)(const double, const double)) {
#pragma omp parallel for schedule(runtime)
  // u(t=0) initial condition
  for (std::ptrdiff_t x = 0; x < mNx; ++x) {
    mU[mNext][x] = pInitU(static_cast<double>(x) * mDx);
  }
  {
    auto time = 0.0;
    auto exact = GetExactSolution(time);
    auto diff = GetDifference(mU[mNext], exact);
    WriteOutput(time, mU[mNext], exact, diff);
  }
  Step();

#pragma omp parallel for schedule(runtime)
  // u_t(t=0) initial conditions
  for (std::ptrdiff_t x = 1; x < mNx - 1; ++x) {
    mU[mNext][x] =
        mU[mCurr][x] + mDt * pInitUT(static_cast<double>(x) * mDx, mC) +
        0.5 * mS2 * (mU[mCurr][x + 1] - 2.0 * mU[mCurr][x] + mU[mCurr][x - 1]);
  }
  mU[mNext][mR] = mU[mCurr][mR] - mS * (mU[mCurr][mR] - mU[mCurr][mR - 1]);
  mU[mNext][0] = mU[mNext][mR];
  {
    auto time = mDt;
    auto exact = GetExactSolution(time);
    auto diff = GetDifference(mU[mNext], exact);
    if (1 % mInterval == 0) {
      WriteOutput(time, mU[mNext], exact, diff);
    }
  }
  Step();
}

void Solver::Solve() {
  for (std::ptrdiff_t t = 2; t < mNt; ++t) {
    SolveOneStep();
    {
      auto time = static_cast<double>(t) * mDt;
      auto exact = GetExactSolution(time);
      auto diff = GetDifference(mU[mNext], exact);
      if (t % mInterval == 0) {
        WriteOutput(time, mU[mNext], exact, diff);
      }
    }
    Step();
  }
}

void Solver::Step() {
  std::swap(mPrev, mCurr);
  std::swap(mCurr, mNext);
}

void Solver::SolveOneStep() {
#pragma omp parallel for schedule(runtime)
  for (std::ptrdiff_t x = 1; x < mNx - 1; ++x) {
    mU[mNext][x] =
        -mU[mPrev][x] + 2.0 * mU[mCurr][x] +
        mS2 * (mU[mCurr][x + 1] - 2.0 * mU[mCurr][x] + mU[mCurr][x - 1]);
  }
  mU[mNext][mR] = mU[mCurr][mR] - mS * (mU[mCurr][mR] - mU[mCurr][mR - 1]);
  mU[mNext][0] = mU[mNext][mR];
}

void Solver::WriteOutput(const double time, const std::vector<double> &rU,
                         const std::vector<double> &rExact,
                         const std::vector<double> &rDiff) {
#pragma omp sections
  {
#pragma omp section
    { WriteToFile(time, rU, mOutputFile); }
#pragma omp section
    { WriteToFile(time, rExact, mExactFile); }
#pragma omp section
    { WriteToFile(time, rDiff, mDiffFile); }
  }
}

std::vector<double> Solver::GetExactSolution(const double time) {
  std::vector<double> result(mNx);
#pragma omp parallel for schedule(runtime)
  for (std::ptrdiff_t x = 0; x < mNx; ++x) {
    result[x] =
        std::sin(2.0 * M_PI * (static_cast<double>(x) * mDx - mC * time));
  }
  return result;
}

std::vector<double> Solver::GetDifference(const std::vector<double> &rArr1,
                                          const std::vector<double> &rArr2) {
  std::vector<double> result(rArr1.size());
#pragma omp parallel for schedule(runtime)
  for (std::size_t i = 0; i < rArr1.size(); ++i) {
    result[i] = rArr1[i] - rArr2[i];
  }
  return result;
}

void Solver::WriteToFile(const double time, const std::vector<double> &rArr,
                         std::ofstream &rOutfile) {
  rOutfile << time;
  for (const auto &val : rArr) {
    rOutfile << "," << val;
  }
  rOutfile << "\n";
}

