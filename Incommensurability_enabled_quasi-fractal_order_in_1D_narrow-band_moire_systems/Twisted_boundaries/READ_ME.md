# Overview

This is a c++ code performs Density Matrix Renormalization Group (DMRG) calculations using the iTensor library (https://itensor.org/).
Includes calculations of the ground state energy, entanglement entropy, fidelity, and structure factor for the generalized Aubry-André-Harper model, using twisted boundary conditions.

## Compilation

The code can simply be compiled using the Makefile provided in this directory via:
  
  ```
  make
  ```

## Running the Code

The code can be ran through:
  ```
  mpiexec -n $ntasks ./dmrg_PD  $id $maxcutoffExp $N $Np $N2 $t $U $V $V2 $parInit $deltaPar $parFinal 
  ```
where ntasks is the number of MPI tasks and the input parameters are described below.  
IMPORTANT: ntasks should be a perfect square number, since the code currently computes
$N'$ different $\phi$ and $N'$ different $k$, so that $ntasks = N'^2$.
This can be changed in the "dmrg_PD.cc" file.

## Input Parameters

The code takes the following input parameters:
1. `id`: Calculation identifier.
2. `maxcutoffExp`: DMRG cutoff exponent. The DMRG cutoff is given by $10^{-\text{maxcutoffExp}}$.
3. `N`: Number of sites.
4. `Np`: Number of particles.
5.  `N2`: Defines $\tau = \frac{N_2}{N}$.
6. `t`: Hopping strength.
7. `U`: Interaction strength.
8. `V`: Aubry-André potential strength.
9. `V2`: Quasiperiodic hopping strength.
10. `parinit`: Initial value of the varying parameter.
11. `deltapar`: Step size for varying the parameter.
12. `parfinal`: Final value of the varying parameter.

Note: To select the varying parameter, one should write "VARY" in the corresponding keyword argument.
Note 2: The $\phi$ and $k$ configurations can be chosen directly in the "dmrg_PD.cc" file


## Output Files

The code generates several output files organized into different folders. 
Below we provide an overview of the output folders and the filenames contained within each folder.

Note: `<filename>` appearing below contains the names and values of all the different parameters, except for the varying parameter.

### Data_Energy Folder

Contains files with ground-state energy values for each $\phi$ and $k$ configuration.

- `<filename>`
  - Data: Each row contains $\phi$ and $k$ values and corresponding energy for each $\phi$ and $k$  configuration.

### Data_Convergence Folder

Contains files with convergence information for each $\phi$ and $k$ configuration.

- `<filename>`
  - Data: Each row contains $\phi$ and $k$ values, maxdim, deltaEnergy, eVariance, deltaS, and NsweepsTot for each $\phi$ and $k$ configuration, where:
    - maxdim: Maximum bond dimension of the MPS
    - deltaEnergy: Difference between the energy of the last two sweeps.
    - eVariance: Energy variance for the last sweep, defined as $\Delta E_{\textrm{var}} = \langle H^2 \rangle - \langle H \rangle^2$.
    - deltaS: Difference between the middle-site entanglement entropies of the last two sweeps.
    - NsweepsTot: Total number of sweeps performed.


### Data_EE Folder

Contains files with entanglement entropy calculations for each $\phi$ and $k$ configuration.

- `<filename>`
  - Data: Each row contains $\phi$ and $k$ values and corresponding middle-site entanglement entropy

### Data_C Folder

Contains files related to the correlation matrix (CMat) and its properties.

- `CMat.<filename>`
  - Data: Matrix elements of the correlation matrix $C_{ij} = \langle c^\dagger_i c_j \rangle$ for each $\phi$ and $k$ configuration.

- `occ.<filename>`
  - Data: $\phi$ value and corresponding single-particle occupancies $\{ n_{\alpha} \}$ for each $\phi$ and $k$ configuration.

- `ni.<filename>`
  - Data: $\phi$ value and corresponding average density $\langle n_i \rangle$ for each $\phi$ and $k$ configuration.

- `S.<filename>`
  - Data: $\phi$ value and corresponding single-particle entanglement entropy defined as $S_{sp} = -\textrm{Tr} (C\log C) = - \sum\limits_{\alpha=1}^N n_{\alpha} \log(n_{\alpha})$ for each $\phi$ and $k$ configuration.

- `IPR.<filename>`
  - Data: $\phi$ value and corresponding IPR of natural orbitals, defined as $IPR_{sp}=\sum\limits_{\alpha=1}^{N} n_{\alpha} \sum\limits_{j=1}^{N} |\phi_{\alpha}(j)|^{4}$, for each $\phi$ and $k$ configuration.

### Data_F Folder

Contains files with fidelity values for each $\phi$ and $k$ configuration.

- `<filename>`
  - Data: Fidelity values for each $\phi$ and $k$ configuration.

### Data_St Folder

Contains files related to the momentum static structure factor ($S(q)$).

- `St.<filename>`
  - Data: Top row contains $q$ values. Bottom row contains $S(q)$ values for each $\phi$ and $k$  configuration. The corresponding value of $\phi$ and $k$ is written at the end of each row.

- `SqPi.<filename>`
  - Data: Each row contains $\phi$ value and corresponding structure factor at $q=\pi$ for each $\phi$ and $k$ configuration.

