/* Program Pairwise PDB, author Ugo Bastolla
   Reads a set of protein sequence and structures.
   Splits pairs of residues into classes (L, [U], |i-j], C^nat)
   Computes a globally optimal selection parameter Lambda by
   minimizing Kullback-Leibler divergence
   Prints scores (d_KL, r^2, RMSE, P(C^nat=1|{a,b}) for each class and globally
   Prints Contact overlap(real, predicted) versus threshold.

*/

int Naa=20; // Do not consider gaps for PDB structures
int NCHAR=1000;

#include "input.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include "coord.h"
#include "protein4.h"
#include "pairwise.h"
#include "Pairwise_aux.h"
#include "optimization.h"
#include "ali_pdb.h"
#include "alignments.h"
#include "REM.h"
//#include "mut_del.h"
//#include "alignments.h"
//#include "gen_code.h"
#include "allocate.h"
#include "read_pdb.h"
//#include "random3.h"           /* Generating random numbers */
//#include "mutation.h"
#include "codes.h"


int ini_print;

/* 
   Give as input a PDB file and either an alignment file or a file with
   mutations. Output is the DeltaG of wild type and mutant sequences.
*/

#define N_CHAR 300          // Max. length of file names


/***********************
    INPUT parameters
************************/
// A: Input files
char name_file[N_CHAR];
char dir_PDB[N_CHAR];
char **file_pdb=NULL, **chain=NULL;
char dir_pdb[300], ext_pdb[10]=".pdb";
char FILE_STR[N_CHAR]="\0", FILE_SEQ[N_CHAR]="\0";
char *file_ali=NULL, *file_mut=NULL;
// B: Thermodynamic parameters
float TEMP=0.5;
//int REM;   // Use 1st (1), 2nd (2) and 3rd (3) moment of misfolded energy 
float **Econt_T=NULL;
// C2: Classes
#define IJ_BIN 40
int N_ij=13, ij_min, ij_max;
int ij_bin[IJ_BIN]={4,5,6,7,8,9,10,15,20,25,30,35,40};
int L_STEP=100, L_min=50, L_max=350;
int N_U=2;
// D: Optimization and weighting
int I_SCORE=0; // Optimize: LIK (0) d_KL (1) or RMSE (2)
//int I_CONT=0;   // Compute the score with 1: C_nat=1 0: C_nat=0,1 2: All 
int I_WEIGHT=0;   // Weight: 0: w=1 1: w=n=Number of pairs 2: w=log(n)


struct pair_class **
Set_pairs_PDB(int *Np_class, int *Nij_class, struct prot_class **prot_class,
	      struct protein *prot_ptr, int N_pdb,
	      int N_ij, int *ij_bin, int N_U,
	      int L_min, int L_max, int L_STEP,
	      float TEMP, int Naa, int I_WEIGHT);
void Print_pairs_PDB(struct prot_class *prot_class, int Np_class,
		     struct pair_class **pairs, int Nij_class,
		     int REM, float TEMP, char *FILE_STR);
void help(char *prog);
int Read_bins(int *bin, char *string);
int Read_pdb_files(char ***file_pdb, char ***chain,
		   char *dir_pdb, char *ext_pdb, char *FILE_IN);
int Get_para(int argc, char **argv,
	     char *FILE_PDB, char *DIR_PDB, float *TEMP, int *REM,
	     char *FILE_STR, char *FILE_SEQ,
	     char **file_ali, char **file_mut,
	     int *ij_bin, int *N_ij, int Nij_max,
	     int *L_min, int *L_max, int *L_STEP,
	     int *N_U, int *I_SCORE, int *I_CONT, int *I_WEIGHT);
extern void Get_name(char *name, char *file_name);
int *Rank(float *X, int n);
int Check_AA(short i_aa, int Naa);

int main(int argc, char **argv){


  /***********************
          DUMMIES
  ************************/
  int i, j;
  char FILE_PDB[300]="\0";

  /***********************
          OUTPUT
  ************************/

  /******************** Input operations   ************************/
  I_CONT=1; REM=2;
  Get_para(argc, argv, FILE_PDB, dir_pdb, &TEMP, &REM, FILE_STR, FILE_SEQ,
	   &file_ali, &file_mut, ij_bin, &N_ij, IJ_BIN, &L_min, &L_max,
	   &L_STEP, &N_U, &I_SCORE, &I_CONT, &I_WEIGHT);
  ij_min=ij_bin[0];
  ij_max=ij_bin[N_ij-1];
  printf("Score to be optimized: ");
  if(I_SCORE==0){printf("Likelihood\n");}
  else if(I_SCORE==1){printf("dKL\n");}
  else if(I_SCORE==2){printf("RMSE\n");}
  printf("Contacts taken into account: ");
  if(I_CONT==0){printf("All\n");}
  else if(I_CONT==1){printf("only native\n");}
  if((I_SCORE==0)&&(I_WEIGHT!=1)){
    printf("WARNING, likelihood score must be used with weights w=n\n");
    I_WEIGHT=1;
  }
  printf("Weight of each class of pairs: ");
  if(I_WEIGHT==0){printf("Equal\n");}
  else if(I_WEIGHT==1){printf("number of pairs\n");}
  else if(I_WEIGHT==2){printf("log(number of pairs)\n");}
  printf("L= %d %d STEP: %d\n", L_min, L_max, L_STEP);
  printf("ij: %d %d Nij= %d\n", ij_min, ij_max, N_ij);
  printf("NU= %d\n", N_U);
  printf("TEMP= %.3f\n", TEMP);

  // Read PDB files
  int N_pdb=0;
  struct protein *prot_ptr=NULL;
  int N_template=0;
  struct protein *template=NULL;
  if(FILE_PDB[0]!='\0'){
    printf("Reading list of PDB files in %s\n", FILE_PDB);
    N_template=Read_pdb_files(&file_pdb, &chain, dir_pdb, ext_pdb, FILE_PDB);
    printf("%d PDB files to read\n", N_template);
    N_template=Read_proteins_PDB(&template, N_template, file_pdb, chain,
				 dir_pdb, ext_pdb, ij_min);
  }
  if((FILE_STR[0]!='\0')&&(FILE_SEQ[0]!='\0')){
    printf("Reading processed proteins in %s and %s\n",
	   FILE_STR, FILE_SEQ);
    N_pdb=Read_processed_proteins(&prot_ptr, FILE_STR, FILE_SEQ, ij_min);
  }else{
    printf("ERROR, no input files specified\n"); exit(8);
  }
  if(N_pdb==0){
    printf("ERROR, no proteins found\n"); exit(8);
  }
  printf("%d proteins read\n", N_pdb);
  
  /* alignment (if any) */
  int PDB=1, i_target=0;
  int n_seq=0, L_ali=0, **ali_seq=NULL, *aligned=NULL;
  char **MSA=NULL, **name_seq=NULL; float **seq_id=NULL;
  short **MSA_aa=NULL, *target=NULL; int L_tar;
  if(file_ali[0]!='\0'){
    float LMIN=0.90; // Minimum fraction of aligned sequence 
    int *selected=NULL;
    MSA=Read_MSA(&n_seq, &L_ali, &name_seq, &selected, file_ali, LMIN);
    if(MSA){
      MSA_aa=Convert_sequences(MSA, n_seq, L_ali);
      ali_seq=Align_pdb(&i_target, &aligned, &seq_id, MSA, n_seq, L_ali,
			file_ali, template, N_template);
      printf("Target sequence: %d\n", i_target);
      L_tar=Remove_gaps_target(MSA_aa, i_target, n_seq, L_ali);
      target=MSA_aa[i_target]; L_max=L_tar;
      PDB=0; // Study alignment instead of PDB
    }
  }

  E_gap=Gap_energy();
  int Indirect=0; // Initialize indirect contacts?
  if(PDB==0){Indirect=1; Naa=21;}
  Econt_T=Energy_over_T(TEMP);

  // Statistics of contacts
  Initialize_REM_2(L_max+1, ij_min, prot_ptr, N_pdb, Indirect);


  int Np_class=0, Nij_class=0;
  struct pair_class **pairs;
  struct prot_class *prot_class;
  if(PDB){ // Indirect contacts are set here!
    pairs=Set_pairs_PDB(&Np_class, &Nij_class, &prot_class, prot_ptr, N_pdb,
			N_ij, ij_bin, N_U, L_min, L_max, L_STEP,
			TEMP, Naa, I_WEIGHT);
  }else{
    float U_ave=Average_energy_over_T(target, L_tar, Econt_T);
    pairs=Set_pairs_ali(&Np_class, &Nij_class, &prot_class,
			U_ave, L_tar, ij_min, Naa);
  }

  // Statistics of indirect contacts
  Indirect_cont_stat(L_max, ij_min, prot_ptr, N_pdb);


  // Optimize Lambda
  if(I_SCORE==0){printf("Optimizing likelihood\n");}
  else if(I_SCORE==1){printf("Optimizing dKL(Q_pred, Q_obs)\n");}
  else if(I_SCORE==2){printf("Optimizing RMSE(Q_pred, Q_obs)\n");}
  else{printf("ERROR, wrong score for optimization\n"); exit(8);}

  int IT_MAX=40, IT_MIN=5; float EPS=0.0001;
  float La_STEP=0.1, La_inf=0.0, La_sup=10, LAMBDA=0.9;
  for(i=0; i<Np_class; i++){
    if(prot_class[i].Np==0)continue;
    float L0=LAMBDA, y0=Score(pairs[i], Naa, Nij_class, I_SCORE, L0);
    float La_opt=L0, y_opt=y0;
    float L1=L0+La_STEP, y1=Score(pairs[i], Naa, Nij_class, I_SCORE, L1);
    if(y1 > y_opt){La_opt=L1; y_opt=y1;}
    float L2=L1+La_STEP, y2=Score(pairs[i], Naa, Nij_class, I_SCORE, L2);
    if(y2 > y_opt){La_opt=L2; y_opt=y2;}
    int iter; float Lambda, y;
    printf("Attempt 0 to optimize Lambda, Lambda=%.3f Score=%.5g\n",
	   La_opt, -y_opt);

    for(iter=0; iter<IT_MAX; iter++){
      Lambda=Find_max_quad(L0, L1, L2, y0, y1, y2, La_inf, La_sup);
      y=Score(pairs[i], Naa, Nij_class, I_SCORE, Lambda);
      if(y > y_opt){La_opt=Lambda; y_opt=y;}
      printf("Attempt %d to optimize Lambda, Lambda=%.3f Score=%.5g\n",
	     iter+1, Lambda, -y);
      if((iter >= IT_MIN)&&(fabs(Lambda-La_opt)<EPS))break; 
      if(Lambda < L0){
	L2=L1; y2=y1; L1=L0; y1=y0; L0=Lambda; y0=y;
      }else if(Lambda < L1){
	L2=L1; y2=y1; L1=Lambda; y1=y;
      }else if(Lambda < L2){
	L0=L1; y0=y1; L1=Lambda; y1=y;
      }else{
	L0=L1; y0=y1; L1=L2; y1=y2; L2=Lambda; y2=y;
      }
    }
    if(iter==IT_MAX)
      printf("WARNING, optimization of Lambda for class %d did not converge\n",
	     i);
    printf("Optimal Lambda: %.3f Score= %.4g\n", La_opt, -y_opt);
    if(fabs(Lambda-La_opt)>EPS)
      y=Score(pairs[i], Naa, Nij_class, I_SCORE, La_opt);
    Prot_stat(prot_class+i, pairs[i], Nij_class, La_opt, Naa);

  }

  // Compute conditional probability of contacts with Bayes theorem
  for(i=0; i<Np_class; i++){
    if(prot_class[i].Np==0)continue;
    for(j=0; j<Nij_class; j++){
      Contact_probability(pairs[i]+j, Naa);
    }
  }

  /*if(file_ali!=NULL){
    Get_name(name_out, file_ali);
    }else if(FILE_PDB[0]!='\0'){
    Get_name(name_out, FILE_PDB);
    }*/
  
  if(1){
    Print_pairs_PDB(prot_class, Np_class, pairs, Nij_class,
		    REM, TEMP, FILE_STR);
  }

  return(0);
}

/************************ Input operations *******************************/

int Get_para(int argc, char **argv,
	     char *FILE_PDB, char *DIR_PDB, float *TEMP, int *REM,
	     char *FILE_STR, char *FILE_SEQ,
	     char **file_ali, char **file_mut,
	     int *ij_bin, int *N_ij, int Nij_max,
	     int *L_min, int *L_max, int *L_STEP, int *N_U,
	     int *I_SCORE, int *I_CONT, int *I_WEIGHT)
{
  int i;
  for(i=1; i<argc; i++){
    if(strncmp(argv[i], "-h", 2)==0)help(argv[0]);
  }
  if(argc < 2){
    printf("ERROR, either input file or options must be specified\n");
    help(argv[0]);
  }

  FILE *file_in=fopen(argv[1], "r");
  if(file_in==NULL){
    printf("WARNING, input file %s does not exist\n", argv[1]);
    return(-1);
  }
  printf("Reading parameters in %s\n", argv[1]);

  char string[1000], READ[100];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    // Protein files
    if(strncmp(string, "FILE_PDB=", 9)==0){
      sscanf(string+9,"%s", FILE_PDB);
   }else if(strncmp(string, "DIR_PDB=", 8)==0){
      sscanf(string+8,"%s", DIR_PDB);
    }else if(strncmp(string, "FILE_STR=", 9)==0){
      sscanf(string+9,"%s", FILE_STR);
    }else if(strncmp(string, "FILE_SEQ=", 9)==0){
      sscanf(string+9,"%s", FILE_SEQ);
    }else if(strncmp(string, "ALI=", 4)==0){
      sscanf(string+4,"%s", READ);
      if(strncmp(READ, "NULL", 4)==0)continue;
      *file_ali=malloc(100*sizeof(char));
      strcpy(*file_ali, READ);
    }else if(strncmp(string, "MUT=", 4)==0){
      sscanf(string+4,"%s", READ);
      if(strncmp(READ, "NULL", 4)==0)continue;
      *file_mut=malloc(100*sizeof(char));
      strcpy(*file_mut, READ);
      // Thermodynamic parameters:
    }else if(strncmp(string, "TEMP=", 5)==0){
      sscanf(string+5,"%f", TEMP);
    }else if(strncmp(string, "REM=", 4)==0){
      sscanf(string+4,"%d", REM);
      if((*REM!=1)&&(*REM!=2)&&(*REM!=3)){
	printf("ERROR, the variable REM can only be set to 1, 2 or 3\n");
	exit(8);
      }
    }else if(strncmp(string, "L_min=", 6)==0){
      sscanf(string+6,"%d", L_min);
    }else if(strncmp(string, "L_max=", 6)==0){
      sscanf(string+6,"%d", L_max);
    }else if(strncmp(string, "L_STEP=", 7)==0){
      sscanf(string+7,"%d", L_STEP);
    }else if(strncmp(string, "N_U=", 4)==0){
      sscanf(string+4,"%d", N_U);
    }else if(strncmp(string, "ij_bin=", 7)==0){
      *N_ij=Read_bins(ij_bin, string+7);
      if(*N_ij>Nij_max){
	printf("ERROR, too many bins of |i-j|, found: %d maximum: %d\n",
	       *N_ij, Nij_max); exit(8);
      }
    }else if(strncmp(string, "Score=", 6)==0){
      sscanf(string+6,"%s", READ);
      if(strncmp(READ, "LIK", 4)==0){*I_SCORE=0;}
      else if(strncmp(READ, "dKL", 3)==0){*I_SCORE=1;}
      else if(strncmp(READ, "RMSE", 4)==0){*I_SCORE=2;}
      else{
	printf("WARNING, unrecognized score %s, using default %d\n",
	       READ, *I_SCORE);
      }
    }else if(strncmp(string, "Cont=", 5)==0){  
     sscanf(string+5,"%d", &i);
     if((i==0)||(i==1)||(i==2)){*I_CONT=i;}
     else{
       printf("WARNING, unrecognized contact type %d, using default %d\n",
	      i, *I_CONT);
     }
    }else if(strncmp(string, "Weight=", 7)==0){  
     sscanf(string+7,"%d", &i);
     if((i==0)||(i==1)||(i==2)){*I_WEIGHT=i;}
     else{
       printf("WARNING, unrecognized weight %d, using default %d\n",
		 i, *I_WEIGHT);
     }
    }else{
      printf("WARNING, uninterpreted line:\n%s", string); 
    }
  }
  fclose(file_in);

  if((FILE_PDB[0]=='\0')&&((FILE_STR[0]=='\0')||(FILE_SEQ[0]=='\0'))){
    printf("ERROR in input file: either list of PDB files (PDB_FILE) or ");
    printf("processed PDBs (FILE_STR, FILE_SEQ) must be specified\n");
    exit(8);
  }

  return(0);
}

int Read_bins(int *bin, char *string){
  char *c=string; int n=0;
  while(*c!='\n'){
    if((*c!=' ')&&(*c!='\t')&&(*c!='\n')){
      sscanf(c, "%d", bin+n); n++;
      while((*c!=' ')&&(*c!='\t')&&(*c!='\n'))c++;
    }
    if(*c=='\n')break; c++;
  }
  return(n);
}

int Read_pdb_files(char ***file_pdb, char ***chain,
		   char *dir_pdb, char *ext_pdb, char *FILE_IN)
{
  int n=0, i;
  FILE *file_in=fopen(FILE_IN, "r");
  char string[1000];
  if(file_in==NULL){
    printf("ERROR, input file %s does not exist\n", FILE_IN);
    return(-1);
  }
  printf("Reading parameters in %s\n", FILE_IN);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue; n++;
  }
  fclose(file_in);
  *file_pdb=malloc(n*sizeof(char *));
  *chain=malloc(n*sizeof(char *));

  file_in=fopen(FILE_IN, "r"); n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "DIR", 3)==0){
      sscanf(string+4, "%s", dir_pdb);
    }else if (strncmp(string, "EXT", 3)==0){
      sscanf(string+4, "%s", ext_pdb);
    }else{
      (*file_pdb)[n]=malloc(80*sizeof(char));
      (*chain)[n]=malloc(10*sizeof(char));
      sscanf(string,"%s%s", (*file_pdb)[n], (*chain)[n]);
      n++;
    }
  }

  // Checks
  for(i=0; i<n; i++){
    if((*file_pdb)[i][0]=='\0'){
      printf("ERROR, PDB file %d not specified\n", i);
    }
  }
  return(n);

}

void help(char *prog){
  printf("USAGE %s <input file> ", prog);
  printf("(ex. Pairwise.in) with default parameters\n");

  printf("\n\nFORMAT of input file:\n");
  printf("#================================================================\n");
  printf("# A) Input files\n");
  printf("PDB_FILE=List_pdb.in or \n");
  printf("FILE_STR=/home/ubastolla/MYPDB/pdbselect25.in  ! Contact matrices\n");
  printf("FILE_SEQ=/home/ubastolla/MYPDB/pdbselect25.seq ! Sequences 1AA per line\n");
  printf("ALI=	family.aln	# FASTA file with MSA (optional)\n");
  //printf("MUT=	prot.mut	# list of mutations (optional)\n");
  printf("#================================================================\n");
  printf("# B) Thermodynamic model\n");
  printf("TEMP=	0.5		# Temperature\n");
  printf("REM=   2		# Use up to 1,2,3 moments of misfold energy\n");
  printf("# C) Site classes\n");
  printf("ij_bin= 4 6 8 10 15 20 25 30 35 40  ! Bins of |i-j|\n");
  printf("L_min=50   # Minimal protein length\n");
  printf("L_max=350  # Maximal protein length\n");
  printf("L_STEP=100 # Size of bin of length\n");
  printf("N_U=3      # Number of bins of average contact energy\n");
  printf("# D) Score computation\n");
  printf("Score=dKL  # Optimized score: LIK, dKL or RMSE\n");
  printf("Cont= 1    # Compute the score: 1: C_nat=1 0: C_nat=0,1 2: All\n");
  printf("Weight=1   # Weight: 0: w=1 1: w=n=Number of pairs 2: w=log(n)\n");
  printf("\n");
  printf("FORMAT of PDB_FILE:\n");
  printf("DIR=/data/ortizg/databases/pdb/ ! Path to pdb files\n");
  printf("EXT=.pdb                        ! Extension of pdb files\n");
  printf("List of pdb codes, one per line\n");
  printf("\n");
  printf("FORMAT of mutation file:\n");
  printf("All mutations that occur together go in the same line\n");
  printf("If more than one chain is present, the chain must be specified");
  printf(" after the mutation.\n");
  printf("If single chain the chain must NOT be specified. Example:\n");
  printf("T103A A T103A B (polychain)\n");
  printf("T103A (single chain)\n");
  printf("FORMAT of deletion:\n");
  printf("DEL T103-W109 A T103A-W109 B (polychain)\n");
  printf("DEL T103-W109 (single chain)\n");

  printf("\n");
  exit(8);
}

int *Rank(float *X, int n){
  // Rank from min (n=0) to max
  int *i_rank=malloc(n*sizeof(int));
  int *ranked=malloc(n*sizeof(int));
  for(int i=0; i<n; i++)ranked[i]=0;
  for(int rank=0; rank<n; rank++){
    int i_min=-1; float X_min;
    for(int i=0; i<n; i++){
      if(ranked[i])continue;
      if((i_min<0)||(X[i]<X_min)){
	i_min=i; X_min=X[i];
      }
    }
    i_rank[rank]=i_min; ranked[i_min]=1;
  }
  free(ranked);
  return(i_rank);
}

int Check_AA(short i_aa, int Naa){
   if((i_aa<0)||(i_aa >=Naa)){
     printf("WARNING, wrong amino acid code %d\n", i_aa);
     return(-1);
   }
   return(0);
 }

struct pair_class **
Set_pairs_PDB(int *Np_class, int *Nij_class, struct prot_class **prot_class,
	      struct protein *prot_ptr, int N_pdb,
	      int N_ij, int *ij_bin, int N_U,
	      int L_min, int L_max, int L_STEP,
	      float TEMP, int Naa, int I_WEIGHT)
{
  // Compute average contact energy per protein
  float U_min=0, U_max=-1000; int i_pdb;
  float *U_ave=malloc(N_pdb*sizeof(float));
  struct protein *prot; int Np=0;
  for(i_pdb=0; i_pdb<N_pdb; i_pdb++){
    prot=prot_ptr+i_pdb;
    if((prot->nres<L_min)||(prot->nres>L_max))continue;
    prot->U_ave=Average_energy_over_T(prot->i_aa, prot->nres, Econt_T);
    U_ave[Np]=prot->U_ave; Np++;
    if(prot->U_ave < U_min)U_min=prot->U_ave;
    if(prot->U_ave > U_max)U_max=prot->U_ave;
  }
  printf("%d proteins\n", Np);
  printf("U_min= %.7f\n", U_min);
  printf("U_max= %.7f\n", U_max);
  int *i_rank=Rank(U_ave, Np); 
  float *U_bin=malloc(N_U*sizeof(float));
  int iU, nU=0, fU=(float)Np/N_U;
  for(iU=0; iU<N_U; iU++){
    nU+=fU; if(nU>=Np)nU=Np-1;
    U_bin[iU]=U_ave[i_rank[nU]];
    printf("Ubin[%d]= %.7f\n", iU, U_bin[iU]);
  }
  free(i_rank);

  /**************************  Allocate  ****************************/
  int N_L, *L_bin;
  struct pair_class **pairs=
    Allocate_pairs(prot_class, Np_class, Nij_class,
		   N_ij, N_U, &N_L, &L_bin, L_min, L_max, L_STEP);
  printf("%d classes of proteins (%d*%d) and ", *Np_class, N_L, N_U);
  printf("%d classes of residue pairs (%d*3)\n", *Nij_class, N_ij);

  // Statistics of pairs
  int Npairs=0;
  int **C_mat=Allocate_mat2_i(L_max, L_max);
  int **C2=Allocate_mat2_i(L_max, L_max);
  int **k_indirect=Allocate_mat2_i(L_max, L_max);
  int **l_indirect=Allocate_mat2_i(L_max, L_max);
  int p_class, ij_class;
  int ij_min=ij_bin[0], ij_max=ij_bin[N_ij-1];
  for(i_pdb=0; i_pdb<N_pdb; i_pdb++){
    prot=prot_ptr+i_pdb;
    if((prot->nres<L_min)||(prot->nres>L_max))continue;
    //if((prot->U_ave<U_min)||(prot->U_ave>U_max))continue;
    p_class=Find_p_class(prot->nres, prot->U_ave, N_L, L_bin, N_U, U_bin);
    struct pair_class *pair_c=pairs[p_class];
    struct prot_class *prot_c=(*prot_class)+p_class;
    prot_c->L+=prot->nres;
    prot_c->U_ave+=prot->U_ave;
    prot_c->Np++;
    Fill_C_mat(C_mat, prot->nres, prot->cont_list, prot->Nc);
    Indirect_contacts(C2, k_indirect, l_indirect, &(prot->indirect_contacts),
		      &prot->Nc_indirect, prot->nres, C_mat);
    for(int i=0; i<prot->nres; i++){
      if(Check_AA(prot->i_aa[i], Naa)<0)continue;
      for(int j=i+ij_min; j<prot->nres; j++){
	if(Check_AA(prot->i_aa[j], Naa)<0)continue;
	Npairs++; int ij=j-i;
	ij_class=Find_ij_class(ij, C_mat[i][j], C2[i][j], N_ij, ij_bin);
	//printf("%d %d\n", p_class, ij_class);
	struct pair_class *pair=pair_c+ij_class;
	pair->ij+=ij;
	pair->P2_obs[prot->i_aa[i]][prot->i_aa[j]]++;
	if(C2[i][j]&&(C_mat[i][j]==0)){ // Indirect contact
	  int k=k_indirect[i][j];
	  if((k<0)||(k>=prot->nres)){
	    printf("ERROR, wrong intermediate residue %d\n", k); exit(8);
	  }
	  pair->P1k[prot->i_aa[k]]++;
	  if((k>i)&&(k<j))pair->k_internal++;
	  pair->lk+=l_indirect[i][j];
	}
      }
    }
  }
  printf("%d pairs stored\n", Npairs);

  int nclass=0;
  for(int i=0; i<(*Np_class); i++){
    struct prot_class *prot=(*prot_class)+i;
    prot->L/=prot->Np;
    prot->U_ave/=prot->Np;
    printf("%d Prot: L= %.0f U= %.5f Np= %d\n",
	   i, prot->L, prot->U_ave, prot->Np);
    if(prot->Np==0)continue;
    double norm=0; int j;
    for(j=0; j<(*Nij_class); j++){
      nclass+=
	Class_statistics(pairs[i]+j, ij_min, ij_max, N_ij, ij_bin, TEMP, Naa);
      if(I_WEIGHT==0){pairs[i][j].w=1;}
      else if(I_WEIGHT==2){pairs[i][j].w=log(pairs[i][j].w);}
      norm+=pairs[i][j].w;
    }
    for(j=0; j<(*Nij_class); j++)pairs[i][j].w/=norm;
  }
  printf("%d classes with data and %d empty\n",
	 nclass, (*Np_class)*(*Nij_class)-nclass);

  return(pairs);
}

void Print_pairs_PDB(struct prot_class *prot_class, int Np_class,
		     struct pair_class **pairs, int Nij_class,
		     int REM, float TEMP, char *FILE_STR)
{
  // Output
  char name_out[300], name_file[300];
  Get_name(name_out, FILE_STR);

  // Write properties of prot_classes: Lambda, d_KL, r^2, RMSE
  sprintf(name_file, "%s_prot_classes_REM%d_T%.2f.dat",
	  name_out, REM, TEMP);
  FILE *file_out=fopen(name_file, "w");
  printf("Writing prot classes in file  %s\n", name_file);
  fprintf(file_out, "#L [U] Lambda Cont lik <d_KL> <RMSE> <r2> n\n");
  int i, j;
  for(i=0; i<Np_class; i++){
    if(prot_class[i].Np==0)continue;
    struct prot_class *prot=prot_class+i;
    for(int c=0; c<3; c++){
      fprintf(file_out, "%.0f\t%.3f\t%.3f\t%d\t%.3e\t%.4f\t%.4f\t%.3f\t%d\n",
	      prot->L, prot->U_ave, prot->Lambda, c, prot->lik[c],
	      prot->d_KL[c], prot->RMSE[c], prot->r2[c], prot->Np);
    }
  }
  fclose(file_out);

  // Write properties of pair_classes: d_KL, r^2, RMSE, P(cont|pairs)
  FILE *file_pair[3];
  for(i=0; i<3; i++){
    sprintf(name_file, "%s_pair_classes_REM%d_T%.2f_Cnat%d.dat",
	    name_out, REM, TEMP, i);
    file_pair[i]=fopen(name_file, "w");
    printf("Writing pair classes with Cnat= %d in file  %s\n", i, name_file);
    fprintf(file_pair[i],
	    "#class L [U] |i-j| Cnat Lambda lik dKL RMSE r2 slope offset P(C|MSA) n\n");
  }
  int k=0;
  for(i=0; i<Np_class; i++){
    struct prot_class *prot=prot_class+i;
    if(prot->Np==0)continue;
    for(j=0; j<Nij_class; j++){
      struct pair_class *pair= pairs[i]+j;
      if(pair->n==0)continue;
      fprintf(file_pair[pair->C_nat],
	      "%d\t%.0f\t%.3f\t%.1f\t%d\t%.3f\t%.2e\t%.4f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d",
	      k, prot->L, prot->U_ave, pair->ij, pair->C_nat, prot->Lambda, 
	      pair->lik, pair->d_KL, pair->RMSE, pair->r2,
	      pair->slope, pair->offset, pair->P_C_nat, pair->n);
      if(pair->C_nat==2)fprintf(file_pair[pair->C_nat],"\t%d",pair->nocont);
      fprintf(file_pair[pair->C_nat],"\n");
      k++;
    }
  }
  for(i=0; i<3; i++)fclose(file_pair[i]);

  // Write Q_obs and Q_pred for each class
  k=0;
  printf("Writing Q_obs and Q_pred in file  %s_Q_k.dat", name_out);
  printf("   0<= k < %d\n", Np_class*Nij_class);
  for(i=0; i<Np_class; i++){
    if(prot_class[i].Np==0)continue;
    for(j=0; j<Nij_class; j++){
      struct pair_class *pair= pairs[i]+j;
      if(pair->n==0)continue;
      sprintf(name_file, "%s_Q_%d.dat", name_out, k);
      file_out=fopen(name_file, "w");
      for(int a=0; a<Naa; a++){
	for(int b=0; b<Naa; b++){
	  fprintf(file_out, "%.3f\t%.3f\t%c%c\n",
		  pair->Q_pred[a][b], pair->Q_obs[a][b],
		  AMIN_CODE[a], AMIN_CODE[b]);
	}
      }
      k++;
      fclose(file_out);
    }
  }


  // Contact overlap
  sprintf(name_file, "%s_contact_overlap_REM%d_T%.2f.dat",
	  name_out, REM, TEMP);
  file_out=fopen(name_file, "w");
  printf("Writing contact overlap in file  %s\n", name_file);
  for(i=0; i<Np_class; i++){
    struct prot_class *prot=prot_class+i;
    struct pair_class *pair=pairs[i], *p;
    for(int c=0; c<3; c++){
      fprintf(file_out, "#%.0f\t%.3f\t%.3f\t%d\t%.4f\t%.4f\t%.3f\t%d\n",
	      prot->L, prot->U_ave, prot->Lambda, c, 
	      prot->d_KL[c]*400, prot->RMSE[c], prot->r2[c], prot->Np);
    }
    float c1=0; for(j=1; j<Nij_class; j+=3)c1+=pair[j].w;
    for(double x=0.1; x<=1.00; x+=0.1){
      float q=0, c2=0; p=pair;
      for(j=0; j<Nij_class; j++){
	if(p->P_C_nat>=x){
	  c2+=p->w; if(p->C_nat==1)q+=p->w;
	}
	p++;
      }
      fprintf(file_out, "%.2f %.3f\n", x, q/sqrt(c1*c2));
    }
    fprintf(file_out, "&\n");
  }
  fclose(file_out);
}
