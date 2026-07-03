int Get_pdb(struct protein *prot, struct residue **res, short **aa_seq,
	    char *file, char *chain, int model);
int Count_models_PDB(char *file_pdb);
int Read_proteins_PDB(struct protein **prot_ptr, int N_prot,
		      char **file_pdb, char **chain,
		      char *dir_pdb, char *ext, int ij_min,
		      int MONOMERIC, int Xray);
int Read_processed_proteins(struct protein **prot_ptr,
			    char *FILE_STR, char *FILE_SEQ,
			    int IJ_MIN, int IJ_MIN_OLD,
			    int MONOMERIC, int Xray, int L_min);
