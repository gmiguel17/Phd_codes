#include <vector>

#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

/*

This file contains functions to compute different physical observables using the
converged DMRG MPS

*/

/* ------------------ FUNCTION TO COMPUTE ENTANGLEMENT ENTROPY
 * -------------------------- */
// Compute entanglement entropy cutting the system between site b and b+1
// (subsystem A: 1 to b; subsytem B: b+1 to N)
auto EntanglementEntropy(auto sitesElec, auto &psi, auto b) {
  //"Gauge" the MPS to site b
  psi.position(b);

  // SVD this wavefunction to get the spectrum
  // of density-matrix eigenvalues
  auto l = leftLinkIndex(psi, b);
  auto s = siteIndex(psi, b);
  auto [U, S, V] = svd(psi(b), {l, s});
  auto u = commonIndex(U, S);

  // Apply von Neumann formula
  // to the squares of the singular values
  Real SvN = 0.;
  for (auto n : range1(dim(u))) {
    auto Sn = elt(S, n, n);
    auto p = sqr(Sn);
    if (p > 1E-12) SvN += -p * log(p);
  }

  return SvN;
}

/* ------------------ FUNCTIONS TO COMPUTE CORRELATION MATRIX
 * -------------------------- */
auto Cij(auto sitesElec, auto &psi, auto i, auto j) {
  auto Adag_i =
      op(sitesElec, "Adag", i);  // automatically returns tensor on site
                                 // i with one prime and one unprimed index
  auto A_j = op(sitesElec, "A", j);

  //'gauge' the MPS to site i
  // any 'position' between i and j, inclusive, would work here
  psi.position(i);

  auto psidag = dag(psi);
  psidag.prime();

  // index linking i to i-1:
  auto li_1 = leftLinkIndex(psi, i);
  auto Cij = prime(psi(i), li_1) * Adag_i * psidag(i);
  for (int k = i + 1; k < j; ++k) {
    Cij *= psi(k);
    Cij *= op(sitesElec, "F", k);  // Jordan-Wigner string
    Cij *= psidag(k);
  }
  // index linking j to j+1:
  auto lj = rightLinkIndex(psi, j);
  Cij *= prime(psi(j), lj);
  Cij *= A_j;
  Cij *= psidag(j);

  // elt (or eltC for complex tensor) extract the elements of a tensor. In this
  // case there is only one element.
  //  see https://itensor.org/docs.cgi?vers=cppv3&page=classes/itensor ("Element
  //  Access Methods")
  auto result = eltC(Cij);  // or eltC(Cij) if expecting complex

  return result;
}

auto Cii(auto sitesElec, auto &psi, auto i) {
  auto N_i = op(sitesElec, "N", i);  // automatically returns tensor on site
                                     // i with one prime and one unprimed index
  //'gauge' the MPS to site i
  psi.position(i);

  auto psidag = dag(psi);
  psidag.prime();

  // index linking i to i-1:
  auto lim = leftLinkIndex(psi, i);
  // index linking i to i+1:
  auto lip = rightLinkIndex(psi, i);
  auto Cii = prime(psi(i), {lim, lip}) * N_i * psidag(i);

  return eltC(Cii);
}

std::vector<Cplx> CijVarj(auto N, auto sitesElec, auto &psi, auto i) {
  std::vector<Cplx> CijElems;
  auto Adag_i =
      op(sitesElec, "Adag", i);  // automatically returns tensor on site
                                 // i with one prime and one unprimed index

  //'gauge' the MPS to site i
  // any 'position' between i and j, inclusive, would work here
  psi.position(i);

  auto psidag = dag(psi);
  psidag.prime();

  // index linking i to i-1:
  auto li_1 = leftLinkIndex(psi, i);
  auto Cij = prime(psi(i), li_1) * Adag_i * psidag(i);
  for (int j = i + 1; j <= N; ++j) {
    if (j > i + 1) {
      Cij *= psi(j - 1);
      Cij *= op(sitesElec, "F", j - 1);  // Jordan-Wigner string
      Cij *= psidag(j - 1);
    }

    auto A_j = op(sitesElec, "A", j);
    auto lj = rightLinkIndex(psi, j);
    // index linking j to j+1:
    auto Cijfinal = Cij * prime(psi(j), lj);
    Cijfinal *= A_j;
    Cijfinal *= psidag(j);
    CijElems.push_back(eltC(Cijfinal));
  }

  return CijElems;
}

auto correlation_matrix(auto sitesElec, auto N, auto &psi) {
  auto ind = Index(N, "ind");
  auto CMat = ITensor(ind, prime(ind));
  for (auto iter = 0; iter < N; iter++) {
    auto cii = Cii(sitesElec, psi, iter + 1);
    CMat.set(ind = iter + 1, prime(ind) = iter + 1, cii);
  }

  for (auto ii = 1; ii <= N; ii++) {
    std::vector<Cplx> Cijelems = CijVarj(N, sitesElec, psi, ii);
    for (auto jj = 0; jj < N - ii; jj++) {
      CMat.set(ind = ii, prime(ind) = ii + jj + 1, Cijelems[jj]);
      CMat.set(ind = ii + jj + 1, prime(ind) = ii, conj(Cijelems[jj]));
    }
  }

  return CMat;
}

auto Cij_MPO(auto sitesElec, auto &psi, auto i, auto j) {
  auto ampo = AutoMPO(sitesElec);
  ampo += 1, "Cdag", i, "C", j;
  auto mpo = toMPO(ampo);
  auto result = innerC(psi, mpo, psi);

  return result;
}

/* ------------------ FUNCTIONS TO COMPUTE DENSITY-DENSITY CORRELATION FUNCTION
 * -------------------------- */
auto ni_nj(auto sitesElec, auto &psi, auto i, auto j) {
  auto N_i = op(sitesElec, "N", i);
  auto N_j = op(sitesElec, "N", j);

  //'gauge' the MPS to site i
  // any 'position' between i and j, inclusive, would work here
  psi.position(i);

  auto psidag = dag(psi);
  psidag.prime();

  // index linking i to i-1:
  auto li_1 = leftLinkIndex(psi, i);
  auto Cij = prime(psi(i), li_1) * (N_i - 1.0 / 2) * psidag(i);
  for (int k = i + 1; k < j; ++k) {
    Cij *= psi(k);
    Cij *= op(sitesElec, "F", k);  // Jordan-Wigner string
    Cij *= psidag(k);
  }
  // index linking j to j+1:
  auto lj = rightLinkIndex(psi, j);
  Cij *= prime(psi(j), lj);
  Cij *= (N_j - 1.0 / 2);
  Cij *= psidag(j);

  // elt (or eltC for complex tensor) extract the elements of a tensor. In this
  // case there is only one element.
  //  see https://itensor.org/docs.cgi?vers=cppv3&page=classes/itensor ("Element
  //  Access Methods")
  auto result = eltC(Cij);  // or eltC(Cij) if expecting complex

  return result;
}

auto ni_ni(auto sitesElec, auto &psi, auto i) {
  auto N_i = op(sitesElec, "N", i);  // automatically returns tensor on site
                                     // i with one prime and one unprimed index

  //'gauge' the MPS to site i
  psi.position(i);

  auto psidag = dag(psi);
  psidag.prime();

  // index linking i to i-1:
  auto lim = leftLinkIndex(psi, i);
  // index linking i to i+1:
  auto lip = rightLinkIndex(psi, i);

  // see http://itensor.org/docs.cgi?vers=julia&page=formulas/mps_onesite_op
  auto IdentityOpi = op(sitesElec, "Id", i);
  auto nini = prime(psi(i), {lim, lip}) * (N_i);

  // Compute <ni>
  auto ni = nini * psidag(i);

  // Compute <ni^2>
  nini.noPrime(inds(N_i)[1]);
  nini *= (N_i)*psidag(i);

  std::vector<Cplx> ni_And_nini;
  ni_And_nini.push_back(eltC(ni));
  ni_And_nini.push_back(eltC(nini));

  return ni_And_nini;
}

std::vector<Cplx> ni_nj_varj(auto N, auto sitesElec, auto &psi, auto i) {
  std::vector<Cplx> CijElems;
  auto N_i = op(sitesElec, "N", i);

  //'gauge' the MPS to site i
  // any 'position' between i and j, inclusive, would work here
  psi.position(i);

  auto psidag = dag(psi);
  psidag.prime();

  // index linking i to i-1:
  // auto IdentityOpi = op(sitesElec,"Id",i);
  auto li_1 = leftLinkIndex(psi, i);
  auto Cij = prime(psi(i), li_1) * (N_i)*psidag(i);
  for (int j = i + 1; j <= N; ++j) {
    if (j > i + 1) {
      Cij *= psi(j - 1);
      Cij *= op(sitesElec, "Id", j - 1);
      Cij *= psidag(j - 1);
    }

    auto N_j = op(sitesElec, "N", j);
    auto lj = rightLinkIndex(psi, j);
    // index linking j to j+1:
    auto Cijfinal = Cij * prime(psi(j), lj);
    // auto IdentityOpj = op(sitesElec,"Id",j);
    Cijfinal *= (N_j);
    Cijfinal *= psidag(j);
    CijElems.push_back(eltC(Cijfinal));
  }

  return CijElems;
}

auto Structure_Factor(auto sitesElec, auto N, auto &psi, auto qlist) {
  auto Npoints = qlist.size();

  // Build (<ninj> - <ni><nj>) matrix
  auto ind = Index(N, "ind");
  auto CMat = ITensor(ind, prime(ind));
  std::vector<Cplx> nivec;
  for (auto iter = 0; iter < N; iter++) {
    auto cii = ni_ni(sitesElec, psi, iter + 1);
    auto diagElem = cii[1] - cii[0] * cii[0];  //<ni^2> - <ni><ni>
    nivec.push_back(cii[0]);
    CMat.set(ind = iter + 1, prime(ind) = iter + 1, diagElem);
  }

  for (auto ii = 1; ii <= N; ii++) {
    std::vector<Cplx> Cijelems = ni_nj_varj(N, sitesElec, psi, ii);

    for (auto jj = 0; jj < N - ii; jj++) {
      auto OffDiagElem =
          Cijelems[jj] - nivec[ii - 1] * nivec[ii + jj];  //<ni nj> - <ni><nj>
      CMat.set(ind = ii, prime(ind) = ii + jj + 1, OffDiagElem);
      CMat.set(ind = ii + jj + 1, prime(ind) = ii,
               OffDiagElem);  // it's the same because  <(ni-1/2)(nj-1/2)>=
                              // <(nj-1/2)(ni-1/2)>
    }
  }

  // Compute structure factors
  auto im = Complex_i;
  std::vector<Cplx> Sqlist;
  for (unsigned int conf = 0; conf < Npoints; conf++) {
    auto Sq = 0.0 + 0.0 * im;
    auto q = qlist[conf];
    for (auto i = 0; i < N; i++) {
      for (auto j = 0; j < N; j++) {
        auto val = eltC(CMat, ind = i + 1, prime(ind) = j + 1);
        val *= (cos(q * (i - j)) + im * sin(q * (i - j)));
        Sq += val;
      }
    }
    Sqlist.push_back(Sq / N / N);
  }

  return Sqlist;
}

auto ni_nj_Mat(auto sitesElec, auto N, auto &psi) {
  // Build (<ninj> - <ni><nj>) matrix
  auto ind = Index(N, "ind");
  auto CMat = ITensor(ind, prime(ind));
  std::vector<Cplx> nivec;
  for (auto iter = 0; iter < N; iter++) {
    auto cii = ni_ni(sitesElec, psi, iter + 1);
    auto diagElem = cii[1] - cii[0] * cii[0];  //<ni^2> - <ni><ni>
    nivec.push_back(cii[0]);
    CMat.set(ind = iter + 1, prime(ind) = iter + 1, diagElem);
  }

  for (auto ii = 1; ii <= N; ii++) {
    std::vector<Cplx> Cijelems = ni_nj_varj(N, sitesElec, psi, ii);

    for (auto jj = 0; jj < N - ii; jj++) {
      auto OffDiagElem =
          Cijelems[jj] - nivec[ii - 1] * nivec[ii + jj];  //<ni nj> - <ni><nj>
      CMat.set(ind = ii, prime(ind) = ii + jj + 1, OffDiagElem);
      CMat.set(ind = ii + jj + 1, prime(ind) = ii,
               OffDiagElem);  // it's the same because  <(ni-1/2)(nj-1/2)>=
                              // <(nj-1/2)(ni-1/2)>
    }
  }

  return CMat;
}

/* ---------------------------------------- */

/* ------------------ FUNCTIONS TO COMPUTE KOHN'S LOCALIZATION TENSOR - AutoMPO
 * Way (just to cross-check) -------------------------- */

auto R_MPO(auto sitesElec, auto &psi, auto N) {
  auto ampo = AutoMPO(sitesElec);
  for (auto i = 1; i <= N; i++) ampo += (i - 1) * 1.0, "Cdag", i, "C", i;
  auto mpo = toMPO(ampo);
  auto result = innerC(psi, mpo, psi);

  return result;
}

auto R2_MPO(auto sitesElec, auto &psi, auto N) {
  auto ampo = AutoMPO(sitesElec);
  for (auto i = 1; i <= N; i++) ampo += (i - 1) * 1.0, "Cdag", i, "C", i;
  auto mpo = toMPO(ampo);

  // Prime MPO A to ensure only one set of site indices are shared
  // see https://itensor.org/docs.cgi?page=classes/mps_mpo_algs
  auto mpo2 = nmultMPO(prime(mpo), mpo);
  auto result = innerC(psi, mpo2, psi);

  return result;
}

auto KohnLocLen_MPO(auto sitesElec, auto &psi, auto N, auto Np) {
  auto psiRpsi = R_MPO(sitesElec, psi, N);
  auto psiR2psi = R2_MPO(sitesElec, psi, N);

  auto loctensor = (psiR2psi - psiRpsi * psiRpsi) / Np;

  if (fabs(imag(loctensor)) >= pow(10, -10.0))
    printf("Problem! Non-real localization length! \n");

  return real(loctensor);
}

/* ---------------------------------------- */

/* ------------------ FUNCTIONS TO COMPUTE KOHN'S LOCALIZATION TENSOR -
 * efficient -------------------------- */

auto xi_xj(auto sitesElec, auto &psi, auto i, auto j) {
  auto N_i = op(sitesElec, "N", i);
  auto N_j = op(sitesElec, "N", j);
  auto x_i = N_i * (i - 1.0);
  auto x_j = N_j * (j - 1.0);

  //'gauge' the MPS to site i
  // any 'position' between i and j, inclusive, would work here
  psi.position(i);

  auto psidag = dag(psi);
  psidag.prime();

  // index linking i to i-1:
  auto li_1 = leftLinkIndex(psi, i);
  auto Cij = prime(psi(i), li_1) * x_i * psidag(i);
  for (int k = i + 1; k < j; ++k) {
    Cij *= psi(k);
    Cij *= op(sitesElec, "F", k);  // Jordan-Wigner string
    Cij *= psidag(k);
  }
  // index linking j to j+1:
  auto lj = rightLinkIndex(psi, j);
  Cij *= prime(psi(j), lj);
  Cij *= x_j;
  Cij *= psidag(j);

  // see https://itensor.org/docs.cgi?vers=cppv3&page=classes/itensor ("Element
  // Access Methods")
  auto result = eltC(Cij);  // or eltC(Cij) if expecting complex

  return result;
}

auto xi_xi(auto sitesElec, auto &psi, auto i) {
  auto N_i = op(sitesElec, "N", i);  // automatically returns tensor on site
                                     // i with one prime and one unprimed index

  auto xi = N_i * (i - 1.0);
  //'gauge' the MPS to site i
  psi.position(i);

  auto psidag = dag(psi);
  psidag.prime();

  // index linking i to i-1:
  auto lim = leftLinkIndex(psi, i);
  // index linking i to i+1:
  auto lip = rightLinkIndex(psi, i);

  // see http://itensor.org/docs.cgi?vers=julia&page=formulas/mps_onesite_op
  auto IdentityOpi = op(sitesElec, "Id", i);
  auto xixi = prime(psi(i), {lim, lip}) * xi;

  // Compute <xi>
  auto xi_av = xixi * psidag(i);

  // Compute <ni^2>
  xixi.noPrime(inds(N_i)[1]);
  xixi *= xi * psidag(i);
  // auto nini = prime(psi(i),{lim,lip})*(N_i*N_i)*psidag(i);

  // println("site i=",i,", nini=",eltC(nini),"\n");

  std::vector<Cplx> xi_And_xixi;
  xi_And_xixi.push_back(eltC(xi_av));
  xi_And_xixi.push_back(eltC(xixi));

  return xi_And_xixi;
}

std::vector<Cplx> xi_xj_varj(auto N, auto sitesElec, auto &psi, auto i) {
  std::vector<Cplx> CijElems;
  auto N_i = op(sitesElec, "N", i);
  auto xi = N_i * (i - 1.0);

  //'gauge' the MPS to site i
  // any 'position' between i and j, inclusive, would work here
  psi.position(i);

  auto psidag = dag(psi);
  psidag.prime();

  // index linking i to i-1:
  // auto IdentityOpi = op(sitesElec,"Id",i);
  auto li_1 = leftLinkIndex(psi, i);
  auto Cij = prime(psi(i), li_1) * (xi)*psidag(i);
  // auto Cij = prime(psi(i),li_1)*N_i*psidag(i) -
  // prime(psi(i),li_1)*IdentityOpi*psidag(i);
  for (int j = i + 1; j <= N; ++j) {
    if (j > i + 1) {
      Cij *= psi(j - 1);
      Cij *= op(sitesElec, "Id", j - 1);
      Cij *= psidag(j - 1);
    }

    auto N_j = op(sitesElec, "N", j);
    auto xj = N_j * (j - 1.0);
    auto lj = rightLinkIndex(psi, j);
    // index linking j to j+1:
    auto Cijfinal = Cij * prime(psi(j), lj);
    // auto IdentityOpj = op(sitesElec,"Id",j);
    Cijfinal *= (xj);
    Cijfinal *= psidag(j);
    CijElems.push_back(eltC(Cijfinal));
  }

  return CijElems;
}

auto KohnLocLen(auto sitesElec, auto N, auto Np, auto &psi) {
  // Build (<xixj> - <xi><xj>) matrix
  auto ind = Index(N, "ind");
  auto CMat = ITensor(ind, prime(ind));
  std::vector<Cplx> nivec;
  for (auto iter = 0; iter < N; iter++) {
    auto cii = xi_xi(sitesElec, psi, iter + 1);
    auto diagElem = cii[1] - cii[0] * cii[0];  //<ni^2> - <ni><ni>
    nivec.push_back(cii[0]);
    CMat.set(ind = iter + 1, prime(ind) = iter + 1, diagElem);
  }

  for (auto ii = 1; ii <= N; ii++) {
    std::vector<Cplx> Cijelems = xi_xj_varj(N, sitesElec, psi, ii);

    for (auto jj = 0; jj < N - ii; jj++) {
      auto OffDiagElem =
          Cijelems[jj] - nivec[ii - 1] * nivec[ii + jj];  //<ni nj> - <ni><nj>
      CMat.set(ind = ii, prime(ind) = ii + jj + 1, OffDiagElem);
      CMat.set(ind = ii + jj + 1, prime(ind) = ii,
               OffDiagElem);  // it's the same because  <(ni-1/2)(nj-1/2)>=
                              // <(nj-1/2)(ni-1/2)>
    }
  }

  // Compute structure factors
  auto im = Complex_i;
  auto kohnloclen = 0.0 + 0.0 * im;
  for (auto i = 0; i < N; i++) {
    for (auto j = 0; j < N; j++) {
      auto val = eltC(CMat, ind = i + 1, prime(ind) = j + 1);
      kohnloclen += val;
    }
  }

  return kohnloclen / Np;
}

/* ------------------ FUNCTIONS TO DUALITY TRANSFORMATION
 * -------------------------- */
// AutoMPS way
auto psiNp_cd_psiNpm1_AUTOMPS(auto sitesElec, auto N, auto &psiNp,
                              auto &psiNpm1) {
  std::vector<Cplx> PHI;
  for (auto i = 1; i <= N; i++) {
    auto sitesElecNoNf = Fermion(N);
    auto ampo = AutoMPO(sitesElecNoNf);
    ampo += "Cdag", i;
    auto mpo = toMPO(ampo);

    auto result = innerC(psiNp, mpo, psiNpm1);
    PHI.push_back(result);
  }

  return PHI;
}

// Efficient way
auto psiNp_cid_psiNpm1(auto sitesElec, auto N, auto &psiNp, auto &psiNpm1,
                       auto i) {
  auto Adag_i =
      op(sitesElec, "Adag", i);  // automatically returns tensor on site
                                 // i with one prime and one unprimed index

  //'gauge' the MPS to site i
  // any 'position' between i and j, inclusive, would work here
  psiNp.position(i);
  psiNpm1.position(i);

  auto psidagNp = dag(psiNp);
  psidagNp.prime();

  auto PHIi = psiNpm1(1);
  if (i == 1) {
    PHIi *= Adag_i;
    PHIi *= psidagNp(1);
  } else {
    PHIi *= op(sitesElec, "Id", 1);
    PHIi *= psidagNp(1);
  }

  for (int k = 2; k < i; ++k) {
    PHIi *= psiNpm1(k);
    PHIi *= op(sitesElec, "Id", k);
    PHIi *= psidagNp(k);
  }

  if (i > 1) {
    PHIi *= psiNpm1(i);
    PHIi *= Adag_i;
    PHIi *= psidagNp(i);
  }
  for (int k = i + 1; k <= N; ++k) {
    PHIi *= psiNpm1(k);
    PHIi *= op(sitesElec, "F", k);  // Jordan-Wigner string
    PHIi *= psidagNp(k);
  }

  // elt (or eltC for complex tensor) extract the elements of a tensor. In this
  // case there is only one element.
  //  see https://itensor.org/docs.cgi?vers=cppv3&page=classes/itensor ("Element
  //  Access Methods")
  auto result = eltC(PHIi);
  println(result);

  return result;
}

auto psiNp_cd_psiNpm1(auto sitesElec, auto N, auto &psiNp, auto &psiNpm1) {
  std::vector<Cplx> PHI;
  for (int i = 1; i <= N; ++i)
    PHI.push_back(psiNp_cid_psiNpm1(sitesElec, N, psiNp, psiNpm1, i));

  return PHI;
}
