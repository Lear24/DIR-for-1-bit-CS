# DIR-for-CS
Distributed Reconstruction via with Invex Relaxation for Compressive Sensing Signals

This repository provides the code for all numerical experiments and empirical experiments in the paper "**Distributed Reconstruction from Compressive Measurements: Nonconvexity and Heterogeneity**".

Due to license restrictions, the raw dataset corresponding to the application in ​**Section 5** of the paper is not provided. However, the experimental code, including data preprocessing, is available.


All necessary custom functions are in the lib folder. When running the code, please add the command:  

`addpath(genpath('.'));`

The simulations contain following  parts:
1. similarity: lambda varies,  given similarity with level of $\pi/3$, $\pi/4$, $\pi/8$.
2. sign-flip probability.
3. Noise variance.
4. Local sample size distribution: The parameter of Dirichlet distribution varies.
5. Varying local sample size with fixed number of nodes.
6. Varying number of nodes with fixed local sample size.
