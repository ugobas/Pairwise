int **Align_pdb(int *target, int **aligned, float ***seq_id,
		char **MSA, int n_seq, int L_ali, char *file_ali,
		struct protein *template, int N_template);
struct pair_class *
Set_pairs_ali(int *Np_class, int *Npairs, struct prot_class **prot_class,
	      int ij_min, int ij_max, int L_ali, int L, struct protein *prot);


int Remove_gaps_target(char **MSA, int i_target, int n_seq, int L_ali,
		       int **ali, struct protein *prots, int N_template,
		       float thr);
void Print_alignments(char *file_ali, int **ali_seq, char **MSA, char **name,
		      int *aligned, struct protein *template, int N_template);

void Copy_MSA(char **MSA, char **MSA_0, int n_seq, int L_ali);
void Convert_sequences(int **MSA_aa, char **MSA, int n_seq, int L_ali);
void Cluster_seqs(int *cluster, char **MSA, int ns, int L_ali, float SI_max);
float *Compute_SI_target(char **MSA, int n_seq, int L_ali, int i_target);
void Weight_seqs(float *w_seq, int *cluster, int ns, float *SI_target, int WSI);
int Select_seqs(char **MSA, int *i_target,
		int n_seq, int L_ali, float *SI_target, float SI_min);
void Count_Pairwise(struct pair_class *pairs, int Npairs,
		    struct prot_class *prot_class, int **MSA_aa,
		    int n_seq, float *w_seq, int L_tar,
		    int ij_min, int Naa); //, int PURGE_REDUNDANCE
void Set_C_nat_ali(struct pair_class *pairs, int Nij_class,
		   struct protein *template, int N_template,
		   int **ali_seq);
