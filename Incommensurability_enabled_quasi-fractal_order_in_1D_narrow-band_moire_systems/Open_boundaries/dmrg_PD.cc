#include <mpi.h>
#include <stdlib.h>

#include <filesystem>
#include <iomanip>  // std::setprecision
#include <vector>

// Do not change the order of the first 4 .h imports!
#include "../Common_code/custom_DMRG.h"  // custom DMRG sweeping procedures
#include "../Common_code/file_handling.h"  // functions to create file names automatically
#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "model_dmrg.h"

using namespace itensor;

int main(int argc, char **argv) {
  // Initialize MPI
  // This must always be called before any other MPI functions
  MPI_Init(&argc, &argv);

  // Get the number of processes in MPI_COMM_WORLD
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of this process in MPI_COMM_WORLD
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // Print out information about MPI_COMM_WORLD
  std::cout << "World Size: " << world_size << "   Rank: " << my_rank
            << std::endl;

  /*  ------------------ Build model and perform DMRG calculation
   * ---------------------- */
  int maximumDim = 10000;  // Very large maxdim --> I will let DMRG use the
                           // largest one needed for the chosen cutoff
  int nsweeps = 7;         // Initial number of sweeps
  auto deltaEnergyCrit = 1E-7;
  auto eVarianceCrit = 1E-6;
  auto deltaSEntCrit = 1E-4;
  auto NsweepsMax = 500;

  // If a given parameter is sent as "VARY", the code assumes that this
  // parameter will be varied Useful for writing to files
  unsigned int inputSize = argc - 1;
  unsigned int exclusion = 100;
  auto flagexclusion = false;
  for (unsigned int iter = 1; iter <= inputSize; iter++) {
    std::string arg = argv[iter];
    if (arg == "VARY") {
      exclusion = iter - 1;
      flagexclusion = true;
    }
  }
  if (flagexclusion == false) {
    println("NO PARAMETER CHOSEN TO VARY!");
  }

  auto id =
      atoi(argv[1]);  // Calculation identifier --> useful to distinguish
                      // between calculations with different internal parameters
  auto maxcutoffExp = atoi(argv[2]);  // DMRG cutoff is 10^(-maxcutoffExp)
  auto maxcutoff = pow(10, -maxcutoffExp);
  auto N = atoi(argv[3]);   // number of sites
  auto Np = atoi(argv[4]);  // number of particles
  auto sitecut =
      int(N / 2.0);  // to compute entanglement entropy
                     //(subsystem A: 1 to sitecut ; subsytem B: sitecut+1 to N)
  auto N2 = atof(argv[5]);  // tau=N2/N
  auto t = atof(argv[6]);   // hopping strength
  auto U = atof(argv[7]);   // interaction strength
  auto V = atof(argv[8]);   // Aubry-André potential strength
  auto V2 = atof(argv[9]);  // Quasiperiodic hopping strength
  unsigned int exclusionparams =
      exclusion -
      2;  // We have id,maxcutoffExp and only then params (SEE BELOW!)
  auto parinit = atof(argv[10]);
  auto deltapar = atof(argv[11]);
  auto parfinal = atof(argv[12]);
  auto Nuc = N;  // number of unit cells in the system (just matters for the
                 // selection of phi configurations below)
  auto dataPrecision =
      12;  // sets the numerical precision of data to write to files

  auto NcTotal = world_size;  // number of phi configurations
  auto phiinit = 3.14159265359 / 4.0;
  auto deltaphi = 2.0 * M_PI / NcTotal / Nuc;
  std::vector<double> phiConfs;
  for (int iphi = 0; iphi < NcTotal; iphi++) {
    phiConfs.push_back(phiinit + iphi * deltaphi);
  }

  // Define sites and sweeps to be used
  auto sitesElec = Fermion(N, {"ConserveNf", true});

  auto sweeps = Sweeps(nsweeps);
  sweeps.maxdim() = 2, 5, 10, 20, 40, 100, maximumDim;
  sweeps.cutoff() = 1E-4, 1E-5, 1E-5, maxcutoff;
  sweeps.noise() = 1E-5, 1E-6, 1E-6, 1E-7, 1E-8, 1E-10;

  // Below we define the configuration for the remaining sweeps that we need to
  // take until the convergence criteria are satisfied
  auto RemainingSweeps = Sweeps(2);
  RemainingSweeps.maxdim() = maximumDim;
  RemainingSweeps.cutoff() = maxcutoff;
  RemainingSweeps.noise() = 1E-10;

  std::cout << parinit << " " << deltapar << " " << parfinal << "\n ";

  MPS psiprev;
  for (auto par = parinit; par <= parfinal + deltapar / 2.0; par += deltapar) {
    // choose configuration according to mpi rank
    auto phi = phiConfs[my_rank];

    std::vector<double> params;
    params.push_back(N);
    params.push_back(Np);
    params.push_back(N2);
    params.push_back(t);
    params.push_back(U);
    params.push_back(V);
    params.push_back(V2);
    params.push_back(phi);
    params[exclusionparams] = par;

    // Definitions for file names
    std::vector<string> VarNames = {"id", "maxCutoff", "N", "Np", "N2",
                                    "t",  "U",         "V", "V2", "phi"};
    std::vector<string> VarValsString = {std::to_string(id),
                                         std::to_string(maxcutoffExp)};
    for (unsigned int iter = 0; iter < params.size(); iter++)
      VarValsString.push_back(std::to_string(params[iter]));
    auto fname = filename(VarNames, VarValsString);

    println("\n **** Variable parameter: ", VarNames[exclusion], "=",
            params[exclusionparams], " **** \n");

    unsigned int phiindexVars = 9;  // index of phi in VarNames
    std::vector<unsigned int> VarExclusions = {
        phiindexVars};  // Exclude positions from 0 to VarNames.size()
    auto fnameExcludepar =
        filenameExcludeParams(VarNames, VarValsString, VarExclusions);

    // Define Hamiltonian
    auto H = model(sitesElec, params);

    // DMRG calculation
    auto [energy, psi, deltaEnergy, eVariance, deltaSEnt, NsweepsTot] =
        DMRGSweepsFullNoInitState(sitesElec, H, params, sweeps, RemainingSweeps,
                                  deltaEnergyCrit, eVarianceCrit, deltaSEntCrit,
                                  NsweepsMax, fname);

    auto maxdim = maxLinkDim(psi);
    println("\n Total QN of Ground State = ",
            totalQN(psi));  // check conservation of total quantum number after
                            // calculation

    // COMPUTE RELEVANT QUANTITIES
    // ///////////////////////////////////////////////////////////////

    /* ----- Compute entanglement entropy ---- */
    auto Sent = EntanglementEntropy(sitesElec, psi, sitecut);

    /* ----- Compute correlation matrix C_ij=<c^\dagger_i c_j> and its
     * eigenvalues and eigenvectors (natural orbitals) ------ */
    auto CMat = correlation_matrix(sitesElec, N, psi);

    // Save CMat in vector<Cplx>
    std::vector<Cplx> CMatVec;
    std::vector<double> niVec;
    auto indicesCMat = inds(CMat);
    auto indRow = indicesCMat[0];
    auto indCol = indicesCMat[1];
    for (int row = 1; row <= N; row++) {
      for (int col = 1; col <= N; col++) {
        auto val = eltC(CMat, indRow = row, indCol = col);
        CMatVec.push_back(val);
        if (row == col) niVec.push_back(val.real());  //<ni> is real anyway!
      }
    }

    // Compute single-particle occupancies and single-particle entropy
    auto [Umat, Dmat] = diagHermitian(CMat);
    auto DmatInd = commonIndex(Umat, Dmat);
    std::vector<double> occupanciesVec;
    for (auto iter = 1; iter <= N; iter++)
      occupanciesVec.push_back(
          elt(Dmat, DmatInd = iter, prime(DmatInd) = iter));

    auto SsinglePart = 0.0;
    for (auto iter = 1; iter <= N; iter++) {
      auto p = elt(Dmat, DmatInd = iter, prime(DmatInd) = iter);
      if (p > 1E-12) SsinglePart += -p * log(p);
    }

    /* Compute IPR of natural orbitals as defined in
     * https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.115.046603 (eq.5) */
    // each column (fixed "UmatColIndex") contains an eigenvector
    auto UmatRowIndex = inds(Umat)[0];
    auto UmatColIndex = inds(Umat)[1];
    auto iprfull = 0.0;
    for (auto evecindex = 1; evecindex <= N; evecindex++) {
      auto ipr = 0.0;
      auto norm = 0.0;
      auto na = elt(Dmat, DmatInd = evecindex, prime(DmatInd) = evecindex);
      for (auto row = 1; row <= N; row++) {
        auto eveci =
            fabs(eltC(Umat, UmatRowIndex = row, UmatColIndex = evecindex));
        ipr += eveci * eveci * eveci * eveci;
        norm += eveci * eveci;
      }
      iprfull += na * ipr;
      if (fabs(norm - 1) > 1.E-8)
        std::cout << "Problem! Eigenvectors of one-particle density matrix not "
                     "normalized! Norm: "
                  << norm << std::endl;
    }
    auto iprfullOverNp = iprfull / Np;

    auto F = 0.0;
    if (par > parinit) {
      F = fabs(innerC(psiprev, psi));
    }

    psiprev = psi;

    /* Compute structure factor */
    std::vector<Cplx> qlist;
    auto Ndivsq = 20;
    auto Npointsq = Ndivsq + 1;
    auto indexPi = 10;  // index of q=pi in qlist
    for (int k = 0; k <= Ndivsq; k++) qlist.push_back(2 * M_PI / Ndivsq * k);
    auto structFactorList = Structure_Factor(sitesElec, N, psi, qlist);
    auto SqPi = structFactorList[indexPi];

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD);  // wait untill all processes reach here!

    // Gathering all data to the root process
    double *phiWorld = NULL;
    std::vector<double> occupanciesVecWorld;
    std::vector<Cplx> CMatVecWorld;
    std::vector<double> nivecWorld;
    std::vector<double> SentlistWorld;
    std::vector<Cplx> structFactorListWorld;
    double *SsinglePartWorld = NULL;
    double *iprfullOverNpWorld = NULL;
    double *FWorld = NULL;
    Cplx *SqPiWorld = NULL;

    double *energyWorld = NULL;
    double *deltaEnergyWorld = NULL;
    double *eVarianceWorld = NULL;
    double *deltaSEntWorld = NULL;
    int *NsweepsTotWorld = NULL;
    int *maxdimWorld = NULL;
    if (my_rank == 0) {
      // Note that in c++ we need to make the casts below, contrary to c!
      phiWorld = (double *)malloc(sizeof(double) * world_size);
      SentlistWorld.resize(world_size);
      occupanciesVecWorld.resize(world_size * N);
      CMatVecWorld.resize(world_size * N * N);
      nivecWorld.resize(world_size * N);
      structFactorListWorld.resize(world_size * Npointsq);
      SsinglePartWorld = (double *)malloc(sizeof(double) * world_size);
      iprfullOverNpWorld = (double *)malloc(sizeof(double) * world_size);
      FWorld = (double *)malloc(sizeof(double) * world_size);
      SqPiWorld = (Cplx *)malloc(sizeof(Cplx) * world_size);

      energyWorld = (double *)malloc(sizeof(double) * world_size);
      deltaEnergyWorld = (double *)malloc(sizeof(double) * world_size);
      eVarianceWorld = (double *)malloc(sizeof(double) * world_size);
      deltaSEntWorld = (double *)malloc(sizeof(double) * world_size);
      NsweepsTotWorld = (int *)malloc(sizeof(int) * world_size);
      maxdimWorld = (int *)malloc(sizeof(int) * world_size);
    }
    MPI_Gather(&phi, 1, MPI_DOUBLE, phiWorld, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&Sent, 1, MPI_DOUBLE, SentlistWorld.data(), 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(occupanciesVec.data(), N, MPI_DOUBLE, occupanciesVecWorld.data(),
               N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(CMatVec.data(), 2 * N * N, MPI_DOUBLE, CMatVecWorld.data(),
               2 * N * N, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);  // the *2 is because CMatVec and CMatVecWorld
                                 // are vectors of complex numbers
    MPI_Gather(niVec.data(), N, MPI_DOUBLE, nivecWorld.data(), N, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(structFactorList.data(), 2 * Npointsq, MPI_DOUBLE,
               structFactorListWorld.data(), 2 * Npointsq, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&SsinglePart, 1, MPI_DOUBLE, SsinglePartWorld, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&iprfullOverNp, 1, MPI_DOUBLE, iprfullOverNpWorld, 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    MPI_Gather(&F, 1, MPI_DOUBLE, FWorld, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&SqPi, 2, MPI_DOUBLE, SqPiWorld, 2, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);  // the 2 is because SqPi is a complex number

    MPI_Gather(&energy, 1, MPI_DOUBLE, energyWorld, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&deltaEnergy, 1, MPI_DOUBLE, deltaEnergyWorld, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&eVariance, 1, MPI_DOUBLE, eVarianceWorld, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&deltaSEnt, 1, MPI_DOUBLE, deltaSEntWorld, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&NsweepsTot, 1, MPI_INT, NsweepsTotWorld, 1, MPI_INT, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&maxdim, 1, MPI_INT, maxdimWorld, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);  // wait till all processes reach here!

    // write to files only from root process after gathering data, to avoid
    // problems
    if (my_rank == 0) {
      std::filesystem::create_directory("Data_Energy");
      std::ofstream outfile;
      outfile.open("Data_Energy/" + fnameExcludepar,
                   std::ios_base::app);  // append instead of overwrite
      for (int ind = 0; ind < world_size; ind++)
        outfile << phiWorld[ind] << " " << std::setprecision(dataPrecision)
                << energyWorld[ind] << "\n";
      outfile.close();

      // Save convergence info
      std::filesystem::create_directory("Data_Convergence");
      outfile.open("Data_Convergence/" + fnameExcludepar,
                   std::ios_base::app);  // append instead of overwrite
      for (int ind = 0; ind < world_size; ind++)
        outfile << phiWorld[ind] << " " << maxdimWorld[ind] << " "
                << std::setprecision(dataPrecision) << deltaEnergyWorld[ind]
                << " " << std::setprecision(dataPrecision)
                << eVarianceWorld[ind] << " "
                << std::setprecision(dataPrecision) << deltaSEntWorld[ind]
                << " " << std::setprecision(dataPrecision)
                << NsweepsTotWorld[ind] << "\n";
      outfile.close();

      // Save entanglement entropy
      std::filesystem::create_directory("Data_EE");
      outfile.open("Data_EE/" + fnameExcludepar,
                   std::ios_base::app);  // append instead of overwrite
      for (int ind = 0; ind < world_size; ind++)
        outfile << phiWorld[ind] << " " << SentlistWorld[ind] << " \n";
      outfile.close();

      // Save CMat -- need to open a different file for each configuration
      std::filesystem::create_directory("Data_C");
      for (int ind1 = 0; ind1 < world_size; ind1++) {
        auto auxVarValsString = VarValsString;
        auxVarValsString[phiindexVars] = std::to_string(phiWorld[ind1]);
        auto fnamevarPhi = filename(VarNames, auxVarValsString);
        outfile.open("Data_C/CMat." +
                     fnamevarPhi);  // I have to use one name per process
        for (int indrow = 0; indrow < N; indrow++) {
          for (int indcol = 0; indcol < N; indcol++) {
            outfile << std::setprecision(dataPrecision)
                    << CMatVecWorld[ind1 * N * N + indrow * N + indcol] << " ";
          }
          outfile << "\n";
        }
        outfile.close();
      }

      // Save occupancies
      outfile.open("Data_C/occ." + fnameExcludepar, std::ios_base::app);
      for (int ind1 = 0; ind1 < world_size; ind1++) {
        outfile << phiWorld[ind1] << " ";
        for (int ind2 = 0; ind2 < N; ind2++) {
          outfile << std::setprecision(dataPrecision)
                  << occupanciesVecWorld[ind1 * N + ind2] << " ";
        }
        outfile << "\n";
      }
      outfile.close();

      // Save <ni>
      outfile.open("Data_C/ni." + fnameExcludepar, std::ios_base::app);
      for (int ind1 = 0; ind1 < world_size; ind1++) {
        outfile << phiWorld[ind1] << " ";
        for (int ind2 = 0; ind2 < N; ind2++) {
          outfile << std::setprecision(dataPrecision)
                  << nivecWorld[ind1 * N + ind2] << " ";
        }
        outfile << "\n";
      }
      outfile.close();

      // Save single particle entanglement entropy
      outfile.open("Data_C/S." + fnameExcludepar,
                   std::ios_base::app);  // append instead of overwrite
      for (int ind = 0; ind < world_size; ind++)
        outfile << phiWorld[ind] << " " << std::setprecision(dataPrecision)
                << SsinglePartWorld[ind] << "\n";
      outfile.close();

      // Save single particle IPR
      outfile.open("Data_C/IPR." + fnameExcludepar,
                   std::ios_base::app);  // append instead of overwrite
      for (int ind = 0; ind < world_size; ind++)
        outfile << phiWorld[ind] << " " << std::setprecision(dataPrecision)
                << iprfullOverNpWorld[ind] << "\n";
      outfile.close();

      // Save Fidelity
      std::filesystem::create_directory("Data_F");
      if (par > parinit) {
        outfile.open("Data_F/" + fnameExcludepar,
                     std::ios_base::app);  // append instead of overwrite
        for (int ind = 0; ind < world_size; ind++)
          outfile << par - deltapar / 2.0 << " " << phiWorld[ind] << " "
                  << std::setprecision(dataPrecision) << FWorld[ind] << "\n";
        outfile.close();
      }

      // Save all structure factors
      std::filesystem::create_directory("Data_St");
      outfile.open("Data_St/" + fnameExcludepar,
                   std::ios_base::app);  // append instead of overwrite
      for (int ind = 0; ind < world_size; ind++) {
        if (ind == 0) {
          for (int indq = 0; indq < Npointsq; indq++)
            outfile << qlist[indq] << " ";
          outfile << "\n";
        }
        for (int indq = 0; indq < Npointsq; indq++)
          outfile << structFactorListWorld[ind * Npointsq + indq] << " ";
        outfile << phiWorld[ind] << " "
                << "\n";
      }
      outfile.close();

      // Save pi structure factor
      outfile.open("Data_St/SqPi." + fnameExcludepar,
                   std::ios_base::app);  // append instead of overwrite
      for (int ind = 0; ind < world_size; ind++) {
        outfile << phiWorld[ind] << " " << std::setprecision(dataPrecision)
                << SqPiWorld[ind] << "\n";
      }
      outfile.close();

      delete[] phiWorld;
      delete[] SsinglePartWorld;
      delete[] iprfullOverNpWorld;
      delete[] FWorld;
      delete[] SqPiWorld;
      delete[] energyWorld;
      delete[] deltaEnergyWorld;
      delete[] eVarianceWorld;
      delete[] deltaSEntWorld;
      delete[] NsweepsTotWorld;
      delete[] maxdimWorld;
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  // Finalize MPI
  MPI_Finalize();
}
