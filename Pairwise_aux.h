//int OPT_MU;     // Optimize mu?
//int I_CONT;  // Compute the score with 1: C_nat=1 0: C_nat=0,1 2: All 


// Energy

float **Energy_over_T(float **E_cont_gap, int Naa, float TEMP);
float **Compute_E_cont_gap(int Naa);
float *Vectorize(float **V2, int Naa);
void Zero_mean(float *V1_not, float *V1,
	       float *P1i1, float *P1i2, //struct pair_class *pair,
	       int Naa);
float Average_energy_over_T(int *i_aa, int L, float T);
struct contact *Indirect_contact_list(int *Nc_indirect,
				      struct contact *cont_list,
				      int Nc, int nres, int ij_min);
void Print_CM(struct protein *prot, FILE *file_out, int *ali);
void Print_Cont_stat(int L, int ij_min);

// Contact statistics
void Compute_Phi(struct pair_class *pair, struct prot_class *prot);
void Compute_ene(struct pair_class *pair, struct prot_class *prot);

void Fill_C_mat(int **C, int nres, struct contact *cont_list, int Nc);
void Indirect_contacts(int **C2, int **k_indirect, int **l_indirect,
		       struct contact **indirect_cont, int *Nc_indirect,
		       int nres, int **C1);
void Indirect_cont_stat(int L_max, int ij_min,
			struct protein *prot_ptr, int N_prot);
float Compute_cont_freq_Nc(int ij, int L);
float Compute_cont_freq(int ij, int L);
float Compute_cont_freq_indir(int ij, int L);
void Initialize_REM_2(int L_max, int ij_min, struct protein *prot,int N_prot);
//void Contact_probability(struct pair_class *pair, int Naa);


// Pairs
void  Initialize_prot_class(struct prot_class *prot);
struct pair_class
**Allocate_pairs_PDB(struct prot_class **prot_class,
		     int *Nprot_class,
		     struct site **sites,
		     int *Nsite_all,
		     int N_C, // Number of sites per protein (contacts)
		     int *Npair_prot, // Number of pairs per protein
		     int N_ij, // number of ranges of i-j
		     int N_U, // NUmber of protein energies
		     int *N_L, // Number of protein sizes
		     int **L_bin, int L_short, int L_long);
void Initialize_pair(struct pair_class *pair, struct prot_class *prot);
void Initialize_site(struct site *site, struct prot_class *prot, int i);
void Empty_pair(struct pair_class *pair);
int Sum_pairs(struct prot_class *prot_c, int Npair_prot,
	      struct protein *protein, int L_seq,
	      int **C_mat, int **C2, int N_ij, int *ij_bin,
	      int ij_min, int N_C, int **label_ij);
struct pair_class *Find_pair(struct pair_class *pair_c, int Npair_prot,
			     int i1, int i2, int C,
			     int ij, int N_ij, int *ij_bin);
void Normalize_sites(struct prot_class *prot, int Nsites_prot, int Npair_prot);

void  Protein_stat(struct prot_class *prot, int Nij_class,
		   int ij_min, int L_max, int I_W);
int Find_p_class(int nres, float U_ave,
		 int N_L, int *L_bin, int N_U, float *U_bin);
int Find_ij_class(int ij, int C_mat, int C2, int N_ij, int *ij_bin,
		  int nc1, int nc2, int N_C, int N_C2);

int Set_index_nc(int nc, int N_C);
int Set_index_nc2(int i1, int i2, int N_C);
int Find_bin_i(int L, int N_L, int *L_bin);
int Find_bin_f(float L, int N_L, float *L_bin);
// Auxiliary
void Copy_vec(float *v1, float *v2, int N);

int Pair_statistics(struct pair_class *pair, int Naa);
float Global_mut_inf(struct pair_class *pair, int Na2);
void Prot_results(struct prot_class *prot_class, int Npair_prot, int L_tar);

void Pairwise_statistics(struct prot_class *prot, int Npair_prot);
int Update_Q_global(struct prot_class *prots, int Np_class,
		     int Npair_prot);
void Compute_Q_global(struct prot_class *prot, int Np, int Npair_prot);

void Symmetrize_mat(double *X);
int Check_nan(double *Q, int n, char *what, int step);
int Check_nan_f(float *Q, int n, char *what, int step);

// Computations
float Score(float Lambda, struct prot_class *prot, int Npair_prot, 
	    int I_SCORE, int comp_all);
/*
void Indirect_Q(float *Qprime, int Naa,
		struct pair_class *pair_ij,
		struct pair_class *pair_ik,
		struct pair_class *pair_jk);
void Indirect_Q_ali(float *Qprime, int Naa, int i_ij,
struct pair_class *pairs);*/

void Find_indirect(struct pair_class *pairs, int Nij_class);

void APC_pairs(struct pair_class *pairs, int Npairs, int L, char *what);
//void Contact_probabilities(struct pair_class *pairs, int Nij_class,
//			   int Naa, int L);
void Contact_probabilities(struct pair_class *pairs, int Nij_class);
