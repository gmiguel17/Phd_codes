#include <vector>

#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

// DMRG calculation starting random initial state
// params[0] and params[1] are ALWAYS assumed to be respectively N and Np
auto dmrgFromModel(auto sitesElec, auto params, auto sweeps) {
  auto H = model(sitesElec, params);

  auto N = int(params[0]);
  auto Np = int(params[1]);

  // Create initial state
  // CDW state until Np is reached
  auto state = InitState(sitesElec);
  auto nocc = 0;
  for (auto i = 1; i <= N;
       i++) {  // i should go from 1 to N (otherwise gives out of range error)
    if (i % 2 == 0 && nocc < Np) {
      state.set(i, "Occ");
      printf("%d ", i);
      nocc++;
    } else
      state.set(i, "Emp");
  }
  // fill the remaining sites
  for (auto i = 1; i <= N;
       i++) {  // i should go from 1 to N (otherwise gives out of range error)
    if (i % 2 != 0 && nocc < Np) {
      state.set(i, "Occ");
      printf("%d ", i);
      nocc++;
    }
  }
  auto psi0 = randomMPS(state);

  auto [energy, psi] =
      dmrg(H, psi0, sweeps,
           {"Quiet", true});  // here I'm using my dmrg function that also
                              // returns the truncation error
  println("\nTotal QN of Ground State = ",
          totalQN(psi));  // check conservation of total quantum number after
                          // calculation
  printfln("Ground state energy = %.20f", energy);

  auto H2 = innerC(
      H, psi, H, psi);  //"innerC" is for complex MPS while "inner" is for real
  auto E = innerC(psi, H, psi);
  auto var = H2 - E * E;
  printfln("Energy check = %.20f", E);
  printfln("Energy variance = %.20f", var);

  return std::make_tuple(energy, psi, var);
}

// DMRG calculation starting with a sent initial state psi0
// params[0] and params[1] are ALWAYS assumed to be respectively N and Np
auto dmrgFromModel(auto sitesElec, auto params, auto sweeps, MPS &psi0) {
  auto H = model(sitesElec, params);

  auto [energy, psi] =
      dmrg(H, psi0, sweeps,
           {"Quiet", true});  // here I'm using my dmrg function that also
                              // returns the truncation error
  println("\nTotal QN of Initial State = ",
          totalQN(psi0));  // check total quantum number
  printfln("Ground state energy = %.20f", energy);

  auto H2 = innerC(
      H, psi, H, psi);  //"innerC" is for complex MPS while "inner" is for real
  auto E = innerC(psi, H, psi);
  auto var = H2 - E * E;
  printfln("Energy check = %.20f", E);
  printfln("Energy variance = %.20f", var);

  return std::make_tuple(energy, psi, var);
}

// Keep sweeping as specified by "RemainingSweeps" untill convergence specified
// by criteria deltaEnergyCrit,deltaSEntCrit and eVarianceCrit is reached If the
// number of sweeps reaches NsweepsMax, the calculation also stops
auto SweepsToConvergence(auto &sitesElec, auto params, auto &RemainingSweeps,
                         auto energy, auto Sent, auto &psi,
                         auto deltaEnergyCrit, auto eVarianceCrit,
                         auto deltaSEntCrit, auto NsweepsMax, auto fname) {
  auto N = int(params[0]);
  auto swstep = 0;
  auto deltaEnergy = deltaEnergyCrit + 0.1;
  auto eVariance = eVarianceCrit + 0.1;
  auto deltaSEnt = deltaSEntCrit + 0.1;
  auto energynew = 0.0;
  auto Sentnew = 0.0;
  while ((deltaEnergy > deltaEnergyCrit || deltaSEnt > deltaSEntCrit ||
          eVariance > eVarianceCrit) &&
         swstep <= NsweepsMax) {
    auto resultnext = dmrgFromModel(sitesElec, params, RemainingSweeps, psi);
    energynew = std::get<0>(resultnext);
    psi = std::get<1>(resultnext);
    eVariance = std::fabs(std::get<2>(resultnext));

    deltaEnergy = std::fabs(energynew - energy);
    energy = energynew;

    Sentnew = EntanglementEntropy(sitesElec, psi, int(N / 2.0));
    deltaSEnt = std::fabs(Sentnew - Sent);
    Sent = Sentnew;

    // println("E=",energy,", S=",Sent,", eVar=",eVariance," \n");
    // println("dE=",deltaEnergy,", dS=",deltaSEnt," \n");

    swstep++;
  }

  return std::make_tuple(energy, psi, deltaEnergy, eVariance, deltaSEnt,
                         swstep);
}

// DMRG calculation starting with random initial MPS (see model_dmrg.cc)
auto DMRGSweepsFullNoInitState(auto &sitesElec, auto params, auto &sweeps,
                               auto &RemainingSweeps, auto deltaEnergyCrit,
                               auto eVarianceCrit, auto deltaSEntCrit,
                               auto NsweepsMax, auto fname) {
  auto N = int(params[0]);
  auto energy = 0.0;
  MPS psi;
  auto deltaEnergy = 0.0;
  auto eVariance = 0.0;
  auto deltaSEnt = 0.0;
  auto NsweepsTot = 0;

  // Initial sweeps
  auto result = dmrgFromModel(sitesElec, params, sweeps);
  energy = std::get<0>(result);
  psi = std::get<1>(result);

  // Remaining sweeps to convergence
  auto Sent = EntanglementEntropy(sitesElec, psi, int(N / 2.0));
  auto final_result = SweepsToConvergence(
      sitesElec, params, RemainingSweeps, energy, Sent, psi, deltaEnergyCrit,
      eVarianceCrit, deltaSEntCrit, NsweepsMax, fname);
  energy = std::get<0>(final_result);
  psi = std::get<1>(final_result);
  deltaEnergy = std::get<2>(final_result);
  eVariance = std::get<3>(final_result);
  deltaSEnt = std::get<4>(final_result);
  NsweepsTot = std::get<5>(final_result);

  return std::make_tuple(energy, psi, deltaEnergy, eVariance, deltaSEnt,
                         NsweepsTot);
}

// This DMRG calculation uses psiold as the initial state for the first sweep
// **IF** loopindex>0 Useful for fidelity calculations, where the MPS in for a
// set of parameters can be used as initial state for the next set of parameters
auto DMRGSweepsFull(auto loopindex, auto &sitesElec, auto params, auto &sweeps,
                    auto &RemainingSweeps, auto psiold, auto deltaEnergyCrit,
                    auto eVarianceCrit, auto deltaSEntCrit, auto NsweepsMax,
                    auto fname) {
  auto N = int(params[0]);
  auto energy = 0.0;
  MPS psi;
  auto deltaEnergy = 0.0;
  auto eVariance = 0.0;
  auto deltaSEnt = 0.0;
  auto NsweepsTot = 0;
  if (loopindex == 0) {
    // Initial sweeps
    auto result = dmrgFromModel(sitesElec, params, sweeps);
    energy = std::get<0>(result);
    psi = std::get<1>(result);

    // Remaining sweeps to convergence
    auto Sent = EntanglementEntropy(sitesElec, psi, int(N / 2.0));
    auto final_result = SweepsToConvergence(
        sitesElec, params, RemainingSweeps, energy, Sent, psi, deltaEnergyCrit,
        eVarianceCrit, deltaSEntCrit, NsweepsMax, fname);
    energy = std::get<0>(final_result);
    psi = std::get<1>(final_result);
    deltaEnergy = std::get<2>(final_result);
    eVariance = std::get<3>(final_result);
    deltaSEnt = std::get<4>(final_result);
    NsweepsTot = std::get<5>(final_result);
  } else {
    // Initial sweeps
    auto result = dmrgFromModel(sitesElec, params, sweeps, psiold);
    energy = std::get<0>(result);
    psi = std::get<1>(result);

    // Remaining sweeps to convergence
    auto Sent = EntanglementEntropy(sitesElec, psi, int(N / 2.0));
    auto final_result = SweepsToConvergence(
        sitesElec, params, RemainingSweeps, energy, Sent, psi, deltaEnergyCrit,
        eVarianceCrit, deltaSEntCrit, NsweepsMax, fname);
    energy = std::get<0>(final_result);
    psi = std::get<1>(final_result);
    deltaEnergy = std::get<2>(final_result);
    eVariance = std::get<3>(final_result);
    deltaSEnt = std::get<4>(final_result);
    NsweepsTot = std::get<5>(final_result);
  }

  return std::make_tuple(energy, psi, deltaEnergy, eVariance, deltaSEnt,
                         NsweepsTot);
}
