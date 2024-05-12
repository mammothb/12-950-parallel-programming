#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>

#include "Solver.h"
#include "Timer.h"
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/parsers.hpp"
#include "boost/program_options/value_semantic.hpp"
#include "boost/program_options/variables_map.hpp"

namespace po = boost::program_options;

/** Initial condition u(x,0) */
double InitU(const double x);
/** Initial condition du(x,t)/dt,t=0 */
double InitUT(const double x, const double c);

int main(int argc, char* argv[]) {
  std::ptrdiff_t nx;        // Mesh points in space
  std::ptrdiff_t nt;        // Mesh points in time
  double c;                 // Constant used to calculate the CFL condition
  std::ptrdiff_t interval;  // Logging interval

  po::options_description description("Program options");
  auto add_opt = description.add_options();
  add_opt("help", "Produce help message");
  add_opt("nx", po::value<std::ptrdiff_t>(&nx)->default_value(100),
          "Mesh points in space");
  add_opt("nt", po::value<std::ptrdiff_t>(&nt)->default_value(200),
          "Mesh points in time");
  add_opt("c", po::value<double>(&c)->default_value(0.1),
          "Constant used to calculate the CFL condition");
  add_opt("interval", po::value<std::ptrdiff_t>(&interval)->default_value(1),
          "Logging interval");

  po::variables_map args;
  po::store(po::parse_command_line(argc, argv, description), args);
  po::notify(args);

  if (args.contains("help")) {
    std::cout << description << std::endl;
    return EXIT_FAILURE;
  }
  if (nx <= 0 || nt <= 0 || interval <= 0) {
    std::cerr << "ERROR: nx, nt, and interval must be positive" << std::endl;
    return EXIT_FAILURE;
  }

  Timer timer;
  Solver solver(nx, nt, c, interval);
  solver.Initialize(&InitU, &InitUT);
  solver.Solve();
  std::cout << "Time: " << timer.Elapsed() << "\n";

  return EXIT_SUCCESS;
}

double InitU(const double x) { return std::sin(2.0 * M_PI * x); }

double InitUT(const double x, const double c) {
  return -2.0 * M_PI * c * std::cos(2.0 * M_PI * x);
}
