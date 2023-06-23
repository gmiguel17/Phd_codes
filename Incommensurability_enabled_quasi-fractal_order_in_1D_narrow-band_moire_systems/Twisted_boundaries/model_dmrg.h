#include <cmath>
#include <vector>

#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

/*

This file contains functions to define the model and carry out a DMRG
calculation. Needs to be changed for different problems

*/

/* ------------------ FUNCTIONS TO BUILD MODEL AND PERFORM DMRG CALCULATIONS
 * -------------------------- */
// params[0] and params[1] are ALWAYS assumed to be respectively N and Np
auto model(auto sitesElec, auto params) {
  // seed for random shuffle that creates initial state
  std::srand(time(NULL));

  auto ii = Complex_i;

  auto N = int(params[0]);
  auto Np = int(params[1]);
  auto N2 = params[2];
  auto tau = (double)N2 / N;
  println("N=", N, ", N2=", N2, " tau=", tau, "\n");
  auto t = params[3];
  auto U = params[4];
  auto V = params[5];
  auto V2 = params[6];
  auto phi = params[7];
  auto k = params[8];  // phase twist

  auto ampo = AutoMPO(sitesElec);
  for (int i = 1; i < N; ++i) {
    ampo += U, "N", i, "N", i + 1;
    ampo += -t * (1.0 + V2 * cos(2 * M_PI * tau * (i - 1 + 1 / 2.0) + phi)),
        "Cdag", i, "C", i + 1;
    ampo += -t * (1.0 + V2 * cos(2 * M_PI * tau * (i - 1 + 1 / 2.0) + phi)),
        "Cdag", i + 1, "C", i;
  }
  for (int i = 1; i <= N; ++i) {
    ampo += V * cos(2 * M_PI * tau * (i - 1) + phi), "N", i;
  }

  // TBC
  ampo += -t * (cos(k) + ii * sin(k)) *
          (1.0 + V2 * cos(2 * M_PI * tau * (N - 1 + 1 / 2.0) + phi)),
      "Cdag", N, "C", 1;
  ampo += -t * (cos(k) - ii * sin(k)) *
          (1.0 + V2 * cos(2 * M_PI * tau * (N - 1 + 1 / 2.0) + phi)),
      "Cdag", 1, "C", N;
  ampo += U, "N", N, "N", 1;

  auto H = toMPO(ampo);

  return H;
}
