# Pairwise

Author Ugo Bastolla ubastolla@cbm.csic.es

## Installation (unix system)

Download Pairwise.zip from  https://github.com/ugobas/Pairwise

Store Pairwise.zip in a new directory

Execute the following:

```sh  
>unzip Pairwise.zip
>make
>cp Pairwise_PDB (PATH) (directory of your executable programs)
```  

## RUN

```sh 
Pairwise_PDB Input_Pairwise.in
``` 

## DESCRIPTION OF THE ALGORITHM

1) Read configuration file Input_Pairwise.in , a set of preprocessed or PDB-file protein strucutures, possibly an alignment

2) Contact statistics:
Compute and store <NC>(L), <C_ij>(l,L), <C_ij NC>(L,l), [U] (for each sequence)

3) Preparation of the data set
3a) Proteins are classified in N_L*N_U classes according to chain length L and average contact energy [U] (hydrophobicity)
3b) Sites of each protein are classified in N_C classes according to their number of contacts
3c) Pairs of sites ij of each protein are classified according to the values of the number of contacts n1, n2 (with n1>=n2 for symmetrizing) and l=|i-j| and native contact C_nat=0 (no contact), C_nat=1 (direct contact) or C_nat=2 (indirect contact). l<3 is excluded and l>l_max is excluded for no contacts. Upper limit for indirect contacts may also be set.

4) For each class of pairs, the program computes the observed pairwise amino-acids distributions P2_ij(a,b) and corresponding single-site P1_i(a)=sum_b P2_ij(a,b) and the observed specific mutual information
Q^obs_ij(a,b)=log(P2_ij(a,b)/(P1_i(a)P1_j(b))).
(it is zero if the two sites i and j are independent, positive if amino acids a and b are positively correlated, otherwise negative)

5) Predicted pairwise distribution based on contact energy U(a,b):
log(Q^pred_ij(a,b))=Lambda*(C^nat_ij-<Cij>+CNC_ij*[U])U(a,b)+eta_i(a)+zeta_j(b)
eta_i(a) and zeta_j(b) are normalization factors, and they are predicted numerically.

The only free parameter is the selection coefficient Lambda, which is determined either by maximizing the likelihood of the observed distribution with respect to the model (Score=LIK in input file), or by minimizing the Kulback-Leibler divergence (Score=dKL), or by minimizg the root mean square error of observed and predicted log(Q) (Score=RME).


## CONFIGURATION FILE

## OUTPUT FILES:
- pdbselect90_pair_classes_<SCORE>_<REM>_<T>_Cnat<N>.dat
Table with the average properties of each set of pair of sites
- pdbselect90_Q_<>.dat
For given pair of sites ij with same number of contacts and different contact distance |i-j| and for each pair of amino acids a,b, table with contact energy and observed and predicted log(Q) for various contact ranges 1i-j| (different columns)
- pdbselect90_prot_classes_<SCORE>_<REM>_<T>.dat
For every one of the N_L*N_U classes of proteins, table with the average properties of each class
- pdbselect90_global_Q.dat
Global pairwise amino acid distribution, averaged over all pairs of sites, and attributed to the mutation process
 

