float Mutate_DG_overT_contfreq(double *Enat, double *E1,
			       double *E1i, double *E2i, 
			       double *E22, double *E23,
			       double *E31, double *E32, double *E323,
			       short *aa_seq, int **C_nat, int L,
			       int res_mut, short aa_new);
float Compute_hydro(short *aa_seq, int L);
float Compute_DG_overT_contfreq(double *Enat, double *E1,
				double *E1i, double *E2i,
				double *E22, double *E23,
				double *E31, double *E32,
				double *E323,
				short *aa_seq, int **C_nat,
				int L, char *file_str,
				float sC0, float sC1,
				float sU1, float T);
float Print_DG_contfreq(double E_nat, double E1, double E22, double E23,
			double E31, double E32, double E323, char *name);
void Test_contfreq(double E1, double E22, double E23, double E31, double E32,
		   double E323, short *aa_seq, int L, char *file_str,
		   char *nameout);
int **Fill_C_nat(int length, short **contact);
int **Get_target(char *file_str, char *name_tar, int *len_tar);
float Compute_Tfreeze(double E1, double E2, double E3, float SC);
void Initialize_REM(float *S_C, float *S_U, int *LEN, int L,
		    char *file_str, float sC0, float sC1, float sU1);
double G_misfold(double E1, double *E2, double *E3,
		 double E22, double E23, double E24,
		 double E32, double E342, double E332, double *E36,
		 float S_C, float *Tf);
float DeltaG(float E_nat, float Gmisf, float S_U);
extern float TEMP, S_C, S_U;
extern int LEN, REM;

//char AA_code[21];
//float Econt[21][21];
//float hydro[20];

//float conf_entropy;
extern double *Cont_freq, NC1, NC2, NC3;
extern double *C221_ij, *C232_i, C242;
extern double *M221_ij, *M232_i, M242;
extern double *C321_ij, *C3X2_ij, *C332_i, C3X3;
extern double *M321_ij, *M3X2_ij, *M332_i, M3X3;
extern double *nc1i, *nc2i, *CNc;
extern double nc2_sum, nc12_sum;
extern int E_36;
extern float **Econt_T, T_ratio;
