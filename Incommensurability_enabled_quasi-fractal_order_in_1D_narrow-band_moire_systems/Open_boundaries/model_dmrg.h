#include <cmath>
#include <vector>

#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

// params[0] and params[1] are ALWAYS assumed to be respectively N and Np
auto model(auto sitesElec, auto params) {
  // seed for random shuffle that creates initial state
  std::srand(time(NULL));

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

  auto H = toMPO(ampo);

  return H;
}
