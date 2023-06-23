# DMRG Code Readme

This code performs a Density Matrix Renormalization Group (DMRG) calculation for a specific model. Below is the documentation for the code.

## Input Parameters

The code takes several input parameters:

1. `id`: Calculation identifier.
2. `maxcutoffExp`: DMRG cutoff exponent. The DMRG cutoff is given by $10^{-\text{maxcutoffExp}}$.
3. `N`: Number of sites.
4. `Np`: Number of particles.
5. $N_2$: Tau parameter, where $\tau = \frac{N_2}{N}$.
6. `t`: Hopping strength.
7. `U`: Interaction strength.
8. `V`: Aubry-André potential strength.
9. `V2`: Quasiperiodic hopping strength.
10. `parinit`: Initial value of the varying parameter.
11. `deltapar`: Step size for varying the parameter.
12. `parfinal`: Final value of the varying parameter.
13. `Nuc`: Number of phi configurations.
14. `dataPrecision`: Numerical precision of data to write to files.

Note: If a given parameter is sent as "VARY", the code assumes that this parameter will be varied. This is useful for writing to files.


## Output Files

The code generates several output files organized into different folders based on the type of data. Here is an overview of the output folders and the files contained within each folder:

### Data_Energy Folder

Contains files with energy values for each $\phi$ configuration.

- `<filename>`
  - Data: $\phi$ value and corresponding energy for each $\phi$ configuration.

### Data_Convergence Folder

Contains files with convergence information for each $\phi$ configuration.

- `<filename>`
  - Data: $\phi$ value, maximum bond dimension, delta energy, energy variance, delta $S_{\text{ent}}$, and total number of sweeps for each $\phi$ configuration.

### Data_EE Folder

Contains files with entanglement entropy values for each $\phi$ configuration.

- `<filename>`
  - Data: $\phi$ value and corresponding entanglement entropy for each $\phi$ configuration.

### Data_C Folder

Contains files related to the correlation matrix (CMat) and its properties.

#### CMat Folder

Contains files with the correlation matrix (CMat) for each $\phi$ configuration.

- `CMat.<filename>`
  - Data: Matrix elements of the correlation matrix $C_{ij} = \langle c^\dagger_i c_j \rangle$ for each $\phi$ configuration.

#### occ Folder

Contains files with single-particle occupancies for each $\phi$ configuration.

- `occ.<filename>`
  - Data: $\phi$ value and corresponding single-particle occupancies for each $\phi$ configuration.

#### ni Folder

Contains files with average particle number $\langle n_i \rangle$ for each $\phi$ configuration.

- `ni.<filename>`
  - Data: $\phi$ value and corresponding average particle number $\langle n_i \rangle$ for each $\phi$ configuration.

#### S Folder

Contains files with single-particle entanglement entropy for each $\phi$ configuration.

- `S.<filename>`
  - Data: $\phi$ value and corresponding single-particle entanglement entropy for each $\phi$ configuration.

#### IPR Folder

Contains files with the inverse participation ratio (IPR) of natural orbitals for each $\phi$ configuration.

- `IPR.<filename>`
  - Data: $\phi$ value and corresponding IPR of natural orbitals for each $\phi$ configuration.

### Data_F Folder

Contains files with fidelity values for each $\phi$ configuration.

- `<filename>`
  - Data: Fidelity values for each $\phi$ configuration.

### Data_St Folder

Contains files related to the structure factor (St) and its properties.

#### St Folder

Contains files with all structure factors ($q$ and corresponding values) for each $\phi$ configuration.

- `St.<filename>`
  - Data: All structure factors ($q$ and corresponding values) for each $\phi$ configuration.

#### SqPi Folder

Contains files with the structure factor at $q=\pi$ for each $\phi$ configuration.

- `SqPi.<filename>`
  - Data: $\phi$ value and corresponding structure factor at $q=\pi$ for each $\phi$ configuration.

Note: `<filename>` is constructed based on the input parameters and varies depending on the specific configuration.