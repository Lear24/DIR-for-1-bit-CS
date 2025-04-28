# DIR-for-CS
Distributed Reconstruction via with Invex Relaxation for Compressive Sensing Signals

This repository provides the code for all numerical experiments and empirical experiments in the paper "**Distributed Reconstruction from Compressive Measurements: Nonconvexity and Heterogeneity**".

All necessary custom functions are in the lib folder. When running the code, please add the command:  

`addpath(genpath('.'));`

The simulations contain following  parts:
1. **test_similarity_pi.m**: similarit; lambda varies,  given similarity with level of $\pi/3$, $\pi/4$, $\pi/8$.
2. **test_sign_filpping_pi.m**: sign-flip probability.
3. **test_noise_variance_pi.m**: Noise variance.
4. **test_local_sample_size_N2400_pi.m**: Local sample size distribution: The parameter of Dirichlet distribution varies.
5. **test_total_sample_vary_m30_pi.m**: Varying local sample size with fixed number of nodes.
6. **test_total_sample_vary_n60_pi.m**: Varying number of nodes with fixed local sample size.



The application on the SEED dataset consists of two parts.
1. **main_final_es_node15_dirichlet_and_nn.m** : Experiment 1 with local measurement quantities in Dirichlet distribution or equal local measurement quantities.
2. **main_final_es_onepatient_dirichlet_and_nn.m** :  Experiment 2 with local measurement quantities in Dirichlet distribution or equal local measurement quantities.

Due to license restrictions, the raw dataset corresponding to the application in ​**Section 5** of the paper is not provided. However, the experimental code, including data preprocessing, is available.
All data processing steps and flowcharts can be found in **Section 5** of the main text and **Appendix C** of the supplementary materials.


Please first run following code after obtain the raw SEED dataset:
`main_empirical_study_patient_is_node.m`
Anyone can reproduce the paper's results by downloading the SEED dataset and running our provided code.
