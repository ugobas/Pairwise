/* Program Pairwise PDB, author Ugo Bastolla
   Reads a set of protein sequence and dstructures.
   Splits pairs of residues into classes (L, [U], |i-j], C^nat)
   Computes a globally optimal selection parameter Lambda by
   minimizing Kullback-Leibler divergence
   Prints scores (dKL, r^2, MSE, P(C^nat=1|{a,b}) for each class and globally
   Prints Contact overlap(real, predicted) versus threshold.

*/
// 0) General
int Naa, Na2;
int PDB;     // Two modes: PDB (statistics of the PDB) or ALI (MSA)
int INDIRECT=1;   // Compute indirect correlations through third residue?
int PRINT_Q=0;   // Print Q matrix?
int Verbose=0;
char dir_pdb[200];
float NTHR; // Smallest number of pairs in a group Na2*FTHR
float FTHR=1.0;

// 1) Classes
// Contact distance
int IJ_MIN=3; // Minimal value of |i-j|
int N_ij=13;
int ij_def[13]={4,5,6,7,8,9,10,15,20,25,30,35,40};
int *ij_bin; // Range of values of |i-j| from input
// Protein length
int L_short=150, L_long=300, L_min=40, L_max=1000;
int N_L=2; // Number of classes of protein length
int N_U=2; // Number of classes of average protein contact energy 
// Number of contacts
int N_C=3; // Classes of number of contacts
int C_low, C_hig; // Few cont: <C_low Many: > C_hig 

// 2) Protein selection
int MONOMERIC=0; // Consider only monomeric proteins?
int Xray=0; //1;      // Consider only Xray structures?

// 3) Misfolding model
int REM=2;     // Use 1st (1), 2nd (2) and 3rd (3) moment of misfolded energy 
float TEMP=1;  // For weighting second moment with REM=2

// 4) Global correlation matrix
int Q_GLOB=1; // Use global Q_matrix (2=ALL_Q)
int IT_GLOB=3; // Number of iterations for optimizing Q_GLOB 
int ALL_Q; // Compute only one Q global for all classes 

// 5) Optimization of Lambda
char Opt_score[40];
int OPT_LAM=1;           // Optimize Lambda?
float La_opt=1.0;        // Value of Lambda if no optimization is performed
int I_SCORE=0;    // Optimize: LIK (0) dKL (1) or MSE (2)
int I_WEIGHT=3;   // Weight: 0: w=1 1: w=n=Numb pairs 2: w=log(n) 3: mut inf
float LAMBDA_MAX=90;  // Maximum value of Lambda

// 6) MSA data
int GAP=1; // Consider gaps in MSA as amino acid 21
int WSI=1;               // Weight sequences by sequence alignment with target
float SI_min=0.25;       // Minimal sequence identity with the target
float SI_max=0.60;       // Cluster and weight seq. if larger SI (Blosum)
float Gap_max=0.25;      // Maximum number of gaps per column admitted
//int PURGE_REDUNDANCE=0;  // Count same pair of a.a. only once? (only for ali)
int N_SI_MIN=3;     // Tested values of SI_min
int N_SI_MAX=3;     // Tested values of SI_max
// WARNING: Number of optimizations of Lambda= N_SI_MIN*N_SI_MAX+IT_GLOB
float Pred_frac=0.333;   // Fraction of contacts to be predicted
int **label_ij; // For alignments, pairs corresponding to a.a. ij i<j

// Common variables
float **E_cont_gap, **E_cont_T, *E_cont_1; // Energy matrix
//float S_DEPEND; // Estimated number of dependent sequences
float *Nc1L, *Nc2L, *Nc3L;

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
#include "REM_new.h"

#include "allocate.h"
#include "read_pdb.h"
//#include "random3.h"           /* Generating random numbers */
//#include "mutation.h"
//#include "mut_del.h"
//#include "alignments.h"
//#include "gen_code.h"
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
char ext_pdb[10]=".pdb";
char FILE_STR[N_CHAR]="\0", FILE_SEQ[N_CHAR]="\0";
char *file_ali=NULL, *file_mut=NULL;
// B: Thermodynamic parameters
//float TEMP=0.5;
//int REM;   // Use 1st (1), 2nd (2) and 3rd (3) moment of misfolded energy 


// Scoring
//float dthr=5.0;

float Score_prot(float *La_opt, struct prot_class *prot,int Npair_prot,
		 int ip);
float Optimize_Lambda(float *La_opt, struct prot_class *prot,
		      int Npair_prot, int i_prot);

struct prot_class *
Set_prots_PDB(int *Np_class, struct site **sites, int N_C,
	      struct pair_class **pairs, int *Npairs,
	      struct protein *prot_ptr, int N_pdb,
	      int N_U, int N_L, int L_min, int L_max, int L_short, int L_long,
	      int N_ij, int *ij_bin,  int Naa);
void Print_pairs(char *name_all, struct prot_class *prot_class, int Np_class,
		 int Npair_prot, int N_C, int PDB, float score_opt);
void Print_average_Q(char *name_all, struct prot_class *prot_class,
		     int Np_class, int Npair_prot);
void Print_log_Q(float *log_Q_pred, float *log_Q_obs, FILE *file_out, int symm);
void Load_score(float *score, struct pair_class *pair, int Npair_prot,int type);

void help(char *prog);
int Read_bins(int *bin, char *string);
int Read_pdb_files(char ***file_pdb, char ***chain,
		   char *dir_pdb, char *ext_pdb, char *FILE_IN);
int Get_para(int argc, char **argv,
	     char *FILE_PDB, char *DIR_PDB, int *MONOMERIC, int *Xray,
	     float *TEMP, int *REM,
	     char *FILE_STR, char *FILE_SEQ,
	     char **file_ali, char **file_mut,
	     int *ij_bin, int *N_ij, int Nij_max,
	     int *L_min, int *L_max, int *L_short, int *L_long,
	     int *N_L, int *N_U, int *N_C, int *Q_GLOB, int *IT_GLOB,
	     int *I_SCORE, int *I_WEIGHT, float *FTHR,
	     float *SI_min, float *SI_max, int *WSI, int *GAP,
	     int *INDIRECT, int *OPT_LAM, int *PRINT_Q);
extern void Get_name(char *name, char *file_name);
int Rank(int *i_rank, float *X, int n);
extern void Sort(int *i_rank, float *X, unsigned long n);

int main(int argc, char **argv){


  /***********************
          DUMMIES
  ************************/
  char FILE_PDB[300]="\0";
  int i;

  /******************** Input operations   ************************/

  ij_bin=malloc(N_ij*sizeof(int));
  for(int i=0; i<N_ij; i++)ij_bin[i]=ij_def[i];
  Get_para(argc, argv, FILE_PDB, dir_pdb, &MONOMERIC, &Xray, &TEMP, &REM,
	   FILE_STR, FILE_SEQ,
	   &file_ali, &file_mut, ij_bin, &N_ij, N_ij, &L_min, &L_max,
	   &L_short, &L_long, &N_L, &N_U, &N_C, &Q_GLOB, &IT_GLOB,
	   &I_SCORE, &I_WEIGHT, &FTHR,
	   &SI_min, &SI_max, &WSI, &GAP, &INDIRECT, &OPT_LAM,
	   &PRINT_Q);
  NTHR=Na2*FTHR;
  if(N_L!=1 && N_L!=2 && N_L!=3){
    printf("ERROR, number of protein length classes N_L=%d "
	   "while it can only be only 1, 2 or 3. Exiting\n", N_L);
    exit(8);
  }
  if(N_U!=1 && N_U!=2 && N_U!=3){
    printf("ERROR, number of protein folding energy classes N_U=%d "
	   "while it can only be only 1, 2 or 3. Exiting\n", N_U);
    exit(8);
  }
  if(N_C!=1 && N_C!=2 && N_C!=3){
    printf("ERROR, number of site classes (number of contacts) N_C=%d "
	   "while it can only be only 1, 2 or 3. Exiting\n", N_C);
    exit(8);
  }
  if(LAMBDA_MAX<=0){
    printf("ERROR, LAMBDA_MAX= %.2g but it must be positive\n", LAMBDA_MAX);
    exit(8);
  }
  if(Q_GLOB==0){IT_GLOB=0;}
  else if(Q_GLOB==1){ALL_Q=1;}
  else if(Q_GLOB==2){ALL_Q=0;}
  else{
    printf("ERROR, not allowed Q_GLOB= %dºn", Q_GLOB); exit(8);
  }

  int ij_max=ij_bin[N_ij-1];
  printf("Score to be optimized: ");
  if(I_SCORE==0){strcpy(Opt_score, "LIK");}
  else if(I_SCORE==1){strcpy(Opt_score, "dKL");}
  else if(I_SCORE==2){strcpy(Opt_score, "MSE");}
  printf("%s\n", Opt_score);
  //printf("Contacts taken into account: ");
  printf("Weight of each class of pairs for opitmizing Lambda: ");
  if(I_WEIGHT==0){printf("Equal\n");}
  else if(I_WEIGHT==1){printf("number of pairs\n");}
  else if(I_WEIGHT==2){printf("log(number of pairs)\n");}
  else if(I_WEIGHT==3){printf("mutual information\n");}
  label_ij=NULL;

  // Read PDB files
  int N_pdb=0;
  struct protein *prot_ptr=NULL;
  if((FILE_STR[0]!='\0')&&(FILE_SEQ[0]!='\0')){
    printf("Reading processed proteins in %s and %s\n",
	   FILE_STR, FILE_SEQ);
    int IJ_MIN_OLD=3;
    N_pdb=Read_processed_proteins(&prot_ptr, FILE_STR, FILE_SEQ,
				  IJ_MIN, IJ_MIN_OLD, MONOMERIC, Xray, L_min);
  }else{
    printf("ERROR, no input files specified\n"); exit(8);
  }
  if(N_pdb==0){
    printf("ERROR, no proteins found\n"); exit(8);
  }
  printf("%d proteins\n", N_pdb);
  struct protein *prot=prot_ptr;
  for(i=0; i<N_pdb; i++){
    prot->indirect_contacts=
      Indirect_contact_list(&prot->Nc_indirect, prot->cont_list,
			    prot->Nc, prot->nres, IJ_MIN);
    prot->min_dist=NULL;
    prot++;
  }
  PDB=1;

  // Read templates, if any
  int N_template=0;
  struct protein *template=NULL;
  if(FILE_PDB[0]!='\0'){
    printf("Reading list of PDB files in %s\n", FILE_PDB);
    N_template=Read_pdb_files(&file_pdb, &chain, dir_pdb, ext_pdb, FILE_PDB);
    printf("%d PDB files to read\n", N_template);
    N_template=Read_proteins_PDB(&template, N_template, file_pdb, chain,
				 dir_pdb, ext_pdb, IJ_MIN, MONOMERIC, Xray);
    struct protein *prot=template;
    for(i=0; i<N_template; i++){
      prot->indirect_contacts=
	Indirect_contact_list(&prot->Nc_indirect, prot->cont_list,
			      prot->Nc, prot->nres, IJ_MIN);
      prot++;
    }
  }

  /* Read alignment (if any) */
  int i_target=-1, rep_pdb=-1, L_pdb=0;
  int n_seq_0=0, L_ali=0, **ali_seq=NULL; // *ali_tar=NULL;
  char **MSA_0=NULL, **name_seq=NULL; float **seq_id=NULL;
  int **MSA_aa=NULL;
  if((file_ali)&&(file_ali[0]!='\0')){
    // ali_seq[i][k]=column aligned to template i site k
    float LMIN=0.50; // Minimum fraction of aligned sequence 
    int *selected=NULL, i;
    MSA_0=Read_MSA(&n_seq_0, &L_ali, &name_seq, &selected, file_ali, LMIN);
    printf("%d sequences read in %s\n", n_seq_0, file_ali);
    if(MSA_0){
      int *aligned=NULL;
      ali_seq=Align_pdb(&i_target, &aligned, &seq_id, MSA_0, n_seq_0, L_ali,
			file_ali, template, N_template);
      if(i_target>=0){    
	printf("Target sequence: %d\n", i_target);
	for(i=0; i<N_template; i++){
	  if((aligned[i]>=0)&&(template[i].nres>L_pdb)){
	    L_pdb=template[i].nres; rep_pdb=i;
	  }
	}
	//if(rep_pdb>=0)ali_tar=ali_seq[rep_pdb];
	printf("Target structure: %d\n", rep_pdb); 
      }else{
	printf("No structure found. Using sequence 0 as target\n");
	i_target=0;
      }
      L_ali=Remove_gaps_target(MSA_0, i_target, n_seq_0, L_ali,
			       ali_seq, template, N_template, Gap_max);
      Print_alignments(file_ali, ali_seq, MSA_0, name_seq, aligned,
		       template, N_template);
      MSA_aa=malloc(n_seq_0*sizeof(int *)); // Numeric codes
      for(int i=0; i<n_seq_0; i++)MSA_aa[i]=malloc(L_ali*sizeof(int));
      Convert_sequences(MSA_aa, MSA_0, n_seq_0, L_ali);
      PDB=0; // Study alignment instead of PDB
    }
  }

  // Print contact-aligned matrices
  int PRINT_CM=0;
  if(PRINT_CM && N_template && ali_seq){
    char nameout[300]; Get_name(nameout, file_ali);
    strcat(nameout, "%s.cm");
    FILE *file_out=fopen(nameout, "w");
    printf("Writing contact matrices in file %s\n", nameout);
    fprintf(file_out, "# MIN contacts. Threshold distance: 4.50 A\n");
    for(i=0; i<N_template; i++){
      Print_CM(template+i, file_out, ali_seq[i]);
    }
    fclose(file_out);
  }
 
/************************   REM computations   ************************/

  Naa=20;  // Do not consider gaps for PDB structures
  if((PDB==0) && GAP)Naa=21; 
  Na2=Naa*Naa;

  // Statistics of contacts
  if(PDB==0)L_max=L_ali;
  if(L_pdb > L_max)L_max=L_pdb;
  Initialize_REM_2(L_max+1, IJ_MIN, prot_ptr, N_pdb);
  // Statistics of indirect contacts
  Indirect_cont_stat(L_max+1, IJ_MIN, prot_ptr, N_pdb);

  if(1){
    int Lp=100; // Print statistics of contacts
    while(Lp<=700){Print_Cont_stat(Lp, IJ_MIN); Lp+=300;}
  }
 
  E_cont_gap=Compute_E_cont_gap(Naa);
  E_cont_T=Energy_over_T(E_cont_gap, Naa, TEMP);
  E_cont_1=Vectorize(E_cont_T, Naa);
  printf("REM= %d TEMP= %.3f", REM, TEMP);
  if(Naa==21)printf(" E_gap=%.3f", E_cont_gap[20][20]);
  printf("\n");

  /******************** Set classes   ************************/
  int Np_class=1, Npair_prot=0, Npairs=0, Nsite_prot=N_C; //, Nsites;
  struct site *sites;
  struct pair_class *pairs;
  struct prot_class *prot_class;
  float score_opt=-1000, score, f_indir=1;

  if(PDB){ 
    /**************** Average over PDB structures *************/
    // Indirect contacts are set here!
    prot_class=Set_prots_PDB(&Np_class, &sites, Nsite_prot, &pairs, &Npairs, 
			     prot_ptr, N_pdb,
			     N_U, N_L, L_min,L_max,L_short,L_long,
			     N_ij, ij_bin, Naa);
    Npair_prot=Npairs/Np_class;
    printf("%d protein classes, %d pairs per protein\n", Np_class, Npair_prot);

    // Npair_prot=3(contacts)*N_C*(N_C+1)/2 (nc)*N_ij
    float Lambda_opt[Np_class];
    for(i=0; i<Np_class; i++){
      struct prot_class *prot=prot_class+i;
      // i=L*N_U+U
      printf("\n### Protein class %d L=%.1f U=%.4f Np=%d\n",
	     i, prot->L/prot->Np, prot->U_ave/prot->U_norm, prot->Np);
      if(prot->Np==0)continue;

      if(ALL_Q==0)Compute_Q_global(prot, 1, Npair_prot);
      Normalize_sites(prot, Nsite_prot, Npair_prot);
      Pairwise_statistics(prot, Npair_prot); // Q_obs computed here
    }
    if(ALL_Q){Compute_Q_global(prot_class, Np_class, Npair_prot);}
    // ALL_Q=1: Same Q_global for all protein classes
    for(i=0; i<Np_class; i++){
      Protein_stat(prot_class+i,Npair_prot,IJ_MIN,L_max,I_WEIGHT);
    }

    int NQ, NP; // Number of Q_glob matrices and proteins per matrix
    if(ALL_Q){NQ=1; NP=Np_class;}else{NQ=Np_class; NP=1;} 
    for(int iq=0; iq<NQ; iq++){
      float Score_max=-1000; int last=0;
      for(int iter=0; iter<=IT_GLOB; iter++){
	printf("### Using Q_global at iteration %d\n", iter);
	float Score_sum=0;
	for(i=0; i<NP; i++){
	  int ip=iq+i;
	  struct prot_class *prot=prot_class+ip;
	  Score_sum+=Optimize_Lambda(&(prot->Lambda), prot, Npair_prot, ip);
	} // end proteins
	printf("### Score at iteration %d: %.4g\n", iter, Score_sum);

	if(iter==0 || Score_sum > Score_max){
	  Score_max=Score_sum;
	  for(i=0; i<NP; i++){
	    int ip=iq+i;
	    struct prot_class *prot=prot_class+ip;
	    Copy_vec(prot->Q_glob_opt, prot->Q_global, Na2);
	    Lambda_opt[ip]=prot->Lambda;
	  }

	}else{ // Score does not improve. Set previous Q_global and leave
	  printf("Score does not improve at iteration %d, "
		 "restoring previous Q_global\n", iter);
	  last=1;
	  for(i=0; i<NP; i++){
	    int ip=iq+i;
	    struct prot_class *prot=prot_class+ip;
	    Copy_vec(prot->Q_global, prot->Q_glob_opt, Na2);
	    //Compute_logarithm(prot->log_Q_glob,prot->Q_global,qmin,lmin,Na2);
	    prot->Lambda=Lambda_opt[ip];
	  }
	}
	if(last || iter==IT_GLOB){
	  printf("Last iteration, computing the scores\n");
	  Score_sum=0;
	  for(i=0; i<NP; i++){
	    int ip=iq+i;
	    struct prot_class *prot=prot_class+ip;
	    float s=Score(prot->Lambda, prot, Npair_prot, I_SCORE, 1);
	    printf("Protein class %d Lambda=%.2f Score=%.6g\n",
		   ip, prot->Lambda, s);
	    Score_sum+=s;
	    Contact_probabilities(prot->pairs, Npair_prot);
	    Prot_results(prot, Npair_prot, N_C);
	  }
	  printf("### Score at last iteration %d: %.4g\n", iter, Score_sum);
	  break;
	} // end last
	int nan=Update_Q_global(prot_class+iq,NP,Npair_prot);
	if(nan){break;}
	// Q_global= sum Q/Q_pred
      } // end iter
      
    } // end iq


    /**************** End optimization of Q_global End PDB *************/
    
  }else{
    /**************** Multiple sequence alignment *************/
    struct protein prot;
    prot.nres=L_max; prot.i_aa=NULL; prot.Nc=Nc1L[L_max]; 
    prot.U_ave=Average_energy_over_T(MSA_aa[i_target], L_ali, TEMP);
    printf("Multiple sequence alignment with %d sequences of length %d\n",
	   n_seq_0, L_ali);
    if(rep_pdb)printf("Representative structure present\n");
    printf("Protein from MSA, L=%d, Nc=%d <U>=%.2f\n",
	   prot.nres, prot.Nc, prot.U_ave);
    pairs=Set_pairs_ali(&Np_class, &Npair_prot, &prot_class,
			IJ_MIN, ij_max, L_ali, L_max, &prot);
    Protein_stat(prot_class,Npair_prot,IJ_MIN,L_max,I_WEIGHT);
    Set_C_nat_ali(pairs, Npair_prot, template, N_template, ali_seq);

 
    // Optimize score with respect to sequence identity
    struct score *scores=malloc(Npair_prot*sizeof(struct score));
    // Parameters
    int n1=N_SI_MAX; float SI_max_ini=SI_max, SI_max_step=0.4;
    int n2=N_SI_MIN; float SI_min_ini=SI_min, SI_min_step=0.1, SI_min_opt=-1;
    int npara=n1*n2;
    // Allocate
    int cluster[n_seq_0];
    float *SI_target_0=Compute_SI_target(MSA_0, n_seq_0, L_ali, i_target);

    char *MSA[n_seq_0]; float SI_target[n_seq_0], w_seq[n_seq_0];
    for(i=0; i<n_seq_0; i++){
      SI_target[i]=SI_target_0[i];
      MSA[i]=malloc(L_ali*sizeof(char));
      Copy_seq(MSA[i], MSA_0[i], L_ali);
    }

    // Loops
    int n_seq=n_seq_0;
    float SI0, SI1, SI2, Sc0, Sc1, Sc2;
    float La_opt_SI=0, f_indir_opt_SI=0, SI_max_opt=-1;
    for(int i1=0; i1<n1; i1++){
      // Exclude sequences k with SI(k,target)<SI_min 
      SI_min=SI_min_ini+i1*SI_min_step;
      n_seq=Select_seqs(MSA, &i_target, n_seq, L_ali, SI_target, SI_min);
      if(n_seq<20)continue;
      float SI_max_tmp=-1, score_tmp=-1; //, f_indir_tmp=0, La_tmp

      for(int i2=0; i2<n2; i2++){
	// Cluster sequences with SI>SI_max
	if(i2<3){
	  SI_max=SI_max_ini+i2*SI_max_step;
	}else{
	  SI_max=Find_max_quad(SI0, SI1, SI2, Sc0, Sc1, Sc2, 0, 1);
	}
	Cluster_seqs(cluster, MSA, n_seq, L_ali, SI_max);
	
	// Weight sequences depending on SI(cl[k],target) if WSI=1
	if(WSI)printf("Weight sequences according to Seq.Id. to target\n");
	Weight_seqs(w_seq, cluster, n_seq, SI_target, WSI);

	Convert_sequences(MSA_aa, MSA, n_seq, L_ali);
	Count_Pairwise(pairs, Npair_prot, prot_class,
		       MSA_aa, n_seq, w_seq, L_ali,
		       IJ_MIN, Naa); //, PURGE_REDUNDANCE
	Compute_Q_global(prot_class, 1, Npair_prot);
	Pairwise_statistics(prot_class, Npair_prot);
	score=Score_prot(&La_opt,prot_class,Npair_prot,0);
	prot_class->Lambda=La_opt;
	if(i2<3){
	  if(i2==0){SI0=SI_max; Sc0=score;}
	  else if(i2==1){SI1=SI_max; Sc1=score;}
	  else if(i2==2){SI2=SI_max; Sc2=score;}
	}else if(SI_max < SI0){
	  SI2=SI1; Sc2=Sc1; SI1=SI0; Sc1=Sc0; SI0=SI_max; Sc0=score;
	}else if(SI_max < SI1){
	    SI2=SI1; Sc2=Sc1; SI1=SI_max; Sc1=score;
	}else if(SI_max < SI2){
	  SI0=SI1; Sc0=Sc1; SI1=SI_max; Sc1=score;
	}else{
	  SI0=SI1; Sc0=Sc1; SI1=SI2; Sc1=Sc2; SI2=SI_max; Sc2=score;
	}
	printf("Sequences: %.0f %.3f < SI < %.3f",
	       pairs[0].n, SI_min, SI_max);
	printf("  Lambda: %.3f Score= %.4g\n", La_opt, score);
	if((i2==0)||(score > score_tmp)){
	  score_tmp=score; f_indir_opt_SI=f_indir;
	  SI_min_opt=SI_min; SI_max_opt=SI_max;
	  if(npara>1)for(int k=0; k<Npair_prot; k++)scores[k]=pairs[k].score;
	}
      } // End loop over SI_max
	if((SI_max_tmp<0)||(score > score_opt)){
	  score_opt=score; La_opt_SI=La_opt; f_indir_opt_SI=f_indir;
	  SI_min_opt=SI_min; SI_max_opt=SI_max;
	  if(npara>1)for(int k=0; k<Npair_prot; k++)scores[k]=pairs[k].score;
	}
    } // End loop over SI_min

    if(npara > 1){
      SI_min=SI_min_opt; SI_max=SI_max_opt;
      La_opt=La_opt_SI; f_indir=f_indir_opt_SI;
      printf("### Optimal parameters: SI_min= %.2f SI_max= %.2f",
	     SI_min, SI_max);
      printf(" Lambda=%.3f f_indir= %.3g Score= %.4g\n",
	     La_opt, f_indir, score_opt);
      Copy_MSA(MSA, MSA_0, n_seq_0, L_ali);
      for(i=0; i<n_seq_0; i++)SI_target[i]=SI_target_0[i];
      // Exclude sequences k with SI(k,target)<SI_min
      n_seq=Select_seqs(MSA, &i_target, n_seq_0, L_ali,SI_target, SI_min);
      Cluster_seqs(cluster, MSA, n_seq_0, L_ali, SI_max);
      // Sequences are weighted depending on SI(cl[k],target) if WSI=1
      if(WSI)printf("Weight sequences according to Seq.Id. to target\n");
      Weight_seqs(w_seq, cluster, n_seq, SI_target, WSI);
      Convert_sequences(MSA_aa, MSA, n_seq, L_ali);
      Count_Pairwise(pairs, Npair_prot, prot_class,
		     MSA_aa, n_seq, w_seq, L_ali, IJ_MIN, Naa);
      Compute_Q_global(prot_class, 1, Npair_prot);
      Pairwise_statistics(prot_class, Npair_prot);
      int Opt_Lam=OPT_LAM; OPT_LAM=0;
      score=Score_prot(&La_opt,prot_class,Npair_prot,0);
      OPT_LAM=Opt_Lam;
    }
    printf("Updating global Q %d times\n", IT_GLOB);
    for(int iter=0; iter<IT_GLOB; iter++){
      Update_Q_global(prot_class, 1, Npair_prot);
      score=Score_prot(&La_opt, prot_class, Npair_prot, 0);
    }
    prot_class->Lambda=La_opt;
    Prot_results(prot_class, Npair_prot, L_ali);
    Contact_probabilities(prot_class->pairs, Npair_prot);
    /***************** End MSA ****************************/
  } 
  
			  
  /***********************  Output  ************************/
  char name_out[100], tmp[20];
   printf("Computations finished, printing output\n");
  if(PDB){Get_name(name_out, FILE_STR);}
  else{Get_name(name_out, file_ali);}
  char name_data[100];
  if(PDB){
    sprintf(name_data, "NL%d_NU%d_NC%d_ij%d", N_L, N_U, N_C, N_ij);
  }else{
    sprintf(name_data, "MSA_");
    if(GAP){strcat(name_data, "GAP");}
    else{strcat(name_data, "noGAP");}
    sprintf(tmp, "_SI%.2f_SI%.2f", SI_min, SI_max);
    strcat(name_data, tmp);
    if(WSI)strcat(name_data, "_WSI");
  }
  char name_model[100];
  sprintf(name_model, "QGLOB%d_REM%d", Q_GLOB, REM);
  if(REM>1){sprintf(tmp, "T%.1f_", TEMP);}else{sprintf(tmp, "_");}
  strcat(name_model, tmp);
  if(OPT_LAM){strcat(name_model, Opt_score);}
  else{sprintf(tmp, "Lamb%.2f", La_opt); strcat(name_model, tmp);}
  
  char name_all[300];
  sprintf(name_all, "%s_", name_out);
  strcat(name_all, name_data);
  strcat(name_all, "_");
  strcat(name_all, name_model);


  Print_pairs(name_all, prot_class, Np_class, Npair_prot, N_C, PDB, score_opt);
  if(PDB==0){
    Print_average_Q(name_all, prot_class, Np_class, Npair_prot);
  }

  return(0);
}

float Score_prot(float *La_opt, struct prot_class *prot, int Npair_prot,
		 int ip)
{
  float score;
  if(OPT_LAM){
    score=Optimize_Lambda(La_opt, prot, Npair_prot, ip);
  }else{
    score=Score(*La_opt, prot, Npair_prot, I_SCORE, 1);
  }
  return(score);
}

/************************ Input operations *******************************/

int Get_para(int argc, char **argv,
	     char *FILE_PDB, char *DIR_PDB, int *MONOMERIC, int *Xray,
	     float *TEMP, int *REM, char *FILE_STR, char *FILE_SEQ,
	     char **file_ali, char **file_mut,
	     int *ij_bin, int *N_ij, int Nij_max,
	     int *L_min, int *L_max, int *L_short, int *L_long,
	     int *N_L, int *N_U, int *N_C, int *Q_GLOB, int *IT_GLOB,
	     int *I_SCORE, int *I_WEIGHT, float *FTHR,
	     float *SI_min, float *SI_max, int *WSI, int *GAP,
	     int *INDIRECT, int *OPT_LAM, int *PRINT_Q)
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

  char string[1000], READ[100]; float xx;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    // Protein files
    if(strncmp(string, "FILE_PDB=", 9)==0){
      sscanf(string+9,"%s", FILE_PDB);
    }else if(strncmp(string, "DIR_PDB=", 8)==0){
      sscanf(string+8,"%s", DIR_PDB);
    }else if(strncmp(string, "MONOMERIC", 9)==0){
      sscanf(string+10,"%d", MONOMERIC);
    }else if(strncmp(string, "XRAY", 4)==0){
      sscanf(string+5,"%d", Xray);
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
    }else if(strncmp(string, "PRINT_Q", 7)==0){
      sscanf(string+8,"%d", PRINT_Q);
      // Thermodynamic parameters:
    }else if(strncmp(string, "TEMP=", 5)==0){
      sscanf(string+5,"%f", TEMP);
    }else if(strncmp(string, "REM=", 4)==0){
      sscanf(string+4,"%d", REM);
      if((*REM<0)&&(*REM>=3)){
	printf("ERROR, the variable REM can only be set to 0, 1, 2 or 3\n");
	exit(8);
      }
    }else if(strncmp(string, "L_min=", 6)==0){
      sscanf(string+6,"%d", L_min);
    }else if(strncmp(string, "L_max=", 6)==0){
      sscanf(string+6,"%d", L_max);
    }else if(strncmp(string, "L_short=", 8)==0){
      sscanf(string+8,"%d", L_short);
      printf("L_short= %d\n", *L_short);
    }else if(strncmp(string, "L_long=", 7)==0){
      sscanf(string+7,"%d", L_long);

    }else if(strncmp(string, "N_L=", 4)==0){
      sscanf(string+4,"%d", N_L);
    }else if(strncmp(string, "N_U=", 4)==0){
      sscanf(string+4,"%d", N_U);
    }else if(strncmp(string, "N_C=", 4)==0){
      sscanf(string+4,"%d", N_C);
    }else if(strncmp(string, "Q_GLOB", 6)==0){
      sscanf(string+7,"%d", Q_GLOB);
    }else if(strncmp(string, "IT_GLOB", 7)==0){
      sscanf(string+8,"%d", IT_GLOB);
    }else if(strncmp(string, "FTHR", 4)==0){
      sscanf(string+5,"%f", &xx);
      if(xx>=0){*FTHR=xx;}
      else{printf("WARNING, FTHR cannot be negative\n");}
      //}else if(strncmp(string, "LAMBDA_MAX", 10)==0){
      //sscanf(string+11,"%f", LAMBDA_MAX);

    }else if(strncmp(string, "ij_bin=", 7)==0){
      *N_ij=Read_bins(ij_bin, string+7);
      if(*N_ij>Nij_max){
	printf("ERROR, too many bins of |i-j|, found: %d maximum: %d\n",
	       *N_ij, Nij_max); exit(8);
      }
    }else if(strncmp(string, "OPT_LAM", 7)==0){
      sscanf(string+8,"%d", OPT_LAM);
    }else if(strncmp(string, "Score=", 6)==0){
      sscanf(string+6,"%s", READ);
      if(strncmp(READ, "LIK", 4)==0){*I_SCORE=0;}
      else if(strncmp(READ, "dKL", 3)==0){*I_SCORE=1;}
      else if(strncmp(READ, "MSE", 4)==0){*I_SCORE=2;}
      else{
	printf("WARNING, unrecognized score %s, using default %d\n",
	       READ, *I_SCORE);
      }
    }else if(strncmp(string, "Weight=", 7)==0){  
     sscanf(string+7,"%d", &i);
     if((i==0)||(i==1)||(i==2)||(i==3)){*I_WEIGHT=i;}
     else{
       printf("WARNING, unrecognized weight %d, using default %d\n",
		 i, *I_WEIGHT);
     }
     printf("I_WEIGHT=%d %s", *I_WEIGHT, string);
    }else if(strncmp(string, "INDIRECT=", 8)==0){  
      sscanf(string+9,"%d", INDIRECT); // Compute indirect contacts
    }else if(strncmp(string, "GAP", 3)==0){
      sscanf(string+4, "%d", GAP);
    }else if(strncmp(string, "Weight_SI=", 9)==0){  
      sscanf(string+10,"%d", WSI);
    }else if(strncmp(string, "SI_min=", 7)==0){  
     sscanf(string+7,"%f", &xx);
     if((xx<0)||(xx>1)){
       printf("WARNING, minimal sequence identity must be between 0 and 1\n");
       printf("Found %.2f, using default %.2f\n", xx, *SI_min);
     }else{
       printf("Setting minimal sequence identity to %.2f\n", xx);
       *SI_min=xx;
     }
    }else if(strncmp(string, "SI_max=", 7)==0){  
     sscanf(string+7,"%f", &xx);
     if((xx<0)||(xx>1)){
       printf("WARNING, maximal sequence identity must be between 0 and 1\n");
       printf("Found %.2f, using default %.2f\n", xx, *SI_max);
     }else{
       printf("Setting maximal sequence identity to %.2f\n", xx);
       *SI_max=xx;
     }
    }else{
      printf("WARNING, uninterpreted line:\n%s", string); 
    }
  }
  fclose(file_in);

  /*
  for(i=1; i<argc; i++){
    if(strcmp(argv[i], "-T", 2)==0){

    }
    }*/


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
    if(*c=='\n'){break;}
    c++;
  }
  if(bin[n-1]<1000){bin[n]=10000; n++;}
  printf("%d bins of ij: ", n);
  for(int i=0; i<n; i++){printf(" %d", bin[i]);}
  printf("\n");
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
    if(string[0]=='#'){continue;} n++;
  }
  fclose(file_in);
  *file_pdb=malloc(n*sizeof(char *));
  *chain=malloc(n*sizeof(char *));

  file_in=fopen(FILE_IN, "r"); n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){continue;}
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
  printf("DIR_PDB= /data/ortizg/databases/pdb/  ! Path to PDB files\n");
  printf("MONOMERIC= 1                          ! Use only monomeric str.?\n");
  printf("XRAY=1                                ! Use only Xray str.?\n");
  printf("PRINT_Q=0                             ! Print Q matrices?\n");
  //printf("MUT=	prot.mut	# list of mutations (optional)\n");
  printf("#================================================================\n");
  printf("# B) Thermodynamic model\n");
  printf("TEMP=	0.5		# Temperature\n");
  printf("REM=   2		# Use up to 1,2,3 moments of misfold energy\n");
  printf("# C) Site classes\n");
  printf("ij_bin= 4 6 8 10 15 20 25 30 35 40  ! Bins of |i-j|\n");
  printf("L_min=30   # Minimal protein length\n");
  printf("L_max=600  # Maximal protein length\n");
  printf("L_short=100 # Short proteins are shorter than this\n");
  printf("L_long=300 # Long proteins are longer than this\n");

  printf("N_U=3      # Number of bins of average contact energy\n");
  printf("# D) Score computation\n");
  printf("OPT_LAM=1  # Optimize Lambda if MSA is given?\n");
  printf("Score=dKL  # Optimized score: LIK, dKL or MSE\n");
  printf("Weight=1   # Weight: 0: w=1 1: w=n=Number of pairs 2: w=log(n)\n");
  printf("INDIRECT=1 # Compute indirect correlations through 1 res? 1=YES\n");
  printf("# E) Alignment options\n");
  printf("GAP=1        # Consider gap as a.a. 20\n");
  printf("SI_min= 0.15 # Minimal sequence identity with the target\n");
  printf("SI_max= 0.60 # Cluster sequences above this identity\n");
  printf("Weight_SI=1  # Weight seqs by Seq.Id. with target 1=YES\n");
  printf("\n");
  printf("FORMAT of PDB_FILE:\n");
  printf("DIR=/ngs/databases/pdb/ ! Path to pdb files\n");
  printf("EXT=.pdb                        ! Extension of pdb files\n");
  printf("List of pdb codes, one per line\n");
  printf("\n");
  /*
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
  */
  printf("\n");
  exit(8);
}

int Rank(int *i_rank, float *X, int n){
  // Rank from min (n=0) to max
  //int *i_rank=malloc(n*sizeof(int));
  int ranked[n], i; for(i=0; i<n; i++)ranked[i]=0;
  for(int rank=0; rank<n; rank++){
    int i_max=-1; float X_max;
    for(i=0; i<n; i++){
      if(ranked[i]){continue;}
      if((i_max<0)||(X[i]>X_max)){
	i_max=i; X_max=X[i];
      }
    }
    i_rank[rank]=i_max; ranked[i_max]=1;
  }
  return(0);
}



struct prot_class *
Set_prots_PDB(int *Np_class, struct site **sites, int N_C,
	      struct pair_class **pairs, int *Npairs,
	      struct protein *prot_ptr, int N_pdb,
	      int N_U, int N_L, int L_min, int L_max, int L_short, int L_long,
	      int N_ij, int *ij_bin,  int Naa)
{

  /**************************  Allocate  ****************************/

  // Allocate proteins
  int Np=N_U*N_L; // Number of protein classes
  struct prot_class *prot_class=malloc(Np*sizeof(struct prot_class));
  *Np_class=Np;

  // Bins of protein length
  int *L_bin=malloc((N_L)*sizeof(int));
  L_bin[0]=L_short;
  if(N_L==3){L_bin[1]=L_long;}
  else if(N_L!=1 && N_L!=2){
    printf("ERROR, N_L=%d while it can only be only 1, 2 or 3. Exiting\n",N_L);
    exit(8);
  }
  L_bin[N_L-1]=100*L_max;

  // Compute average mean energy per protein and contacts per residue
  int L_max_act=0; // Largest found protein
  double U_mean=0, nc_mean=0;
  float U_min=1000, U_max=-1000; 
  int nc_min=1000, nc_max=-1;
  Np=0;
  for(int i_pdb=0; i_pdb<N_pdb; i_pdb++){
    struct protein *prot=prot_ptr+i_pdb;
    if((prot->nres<L_min)||(prot->nres>L_max)){continue;}
    Np++;
    if(prot->nres>L_max_act){L_max_act=prot->nres;}
    prot->U_ave=Average_energy_over_T(prot->i_aa, prot->nres, TEMP);
    if(prot->U_ave < U_min)U_min=prot->U_ave;
    if(prot->U_ave > U_max)U_max=prot->U_ave;
    float sum_nc=0;
    for(int k=0; k<prot->nres; k++){
      int nc=prot->n_cont[k]; sum_nc+=nc;
      if(nc<nc_min){nc_min=nc;}
      if(nc>nc_max){nc_max=nc;}
    }
    nc_mean+=sum_nc/prot->nres;
  }
  nc_mean/=Np; U_mean/=Np;
  printf("%d proteins\n", Np);

  // Bins of folding free energy
  float *U_bin=malloc(N_U*sizeof(float));
  float U_low, U_high;
  if(N_U==2){
    U_low=U_mean; U_high=U_max;
    U_bin[0]=U_low; U_bin[1]=U_high;
  }else if(N_U==3){
    float U_range=(U_max-U_min)/3;
    U_low=U_min+U_range;
    U_high=U_max-U_range;
    U_bin[0]=U_low; U_bin[1]=U_high; U_bin[2]=2*U_max;
  }else if(N_U==1){
    U_low=U_max; U_high=2*U_max;
    U_bin[0]=U_max;
  }else{
    printf("ERROR, wrong number of bins for energy\n"); exit(8);
  }
  printf("energy per protein: min= %.3g mean= %.3g max= %.3g\n",
	 U_min, U_mean, U_max);
  printf("%d energy bins: < %.3g > %.3g\n", N_U, U_low, U_high);

  // Set sites
  int Nsite_all=Np*N_C;
  *sites=malloc(Nsite_all*sizeof(struct site));
  if(N_C==3){
    float nc_third=(nc_max-nc_min)/3;
    C_low=nc_min+nc_third;
    C_hig=nc_max-nc_third;
  }else if(N_C==2){
    C_low=nc_mean;
    C_hig=nc_max;
  }else if(N_C==1){
    C_low=nc_max;
    C_hig=2*nc_max;
  }
  printf("%d  classes of number of contacts: mean %.1f low < %d high > %d\n",
	 N_C, nc_mean, C_low, C_hig);

  // Allocate pairs
  int Npair_sites=N_C*(N_C+1)/2; // (number of pairs of sites)
  int Npair_prot=3*Npair_sites*N_ij; // NUmber of pairs per protein
  *Npairs=*Np_class*Npair_prot;
  *pairs=malloc((*Npairs)*sizeof(struct pair_class));

  printf("%d classes of proteins (%d*%d) and ", *Np_class, N_L, N_U);
  printf("%d classes of residue pairs "
	 "(3 cont_types *%d site_pairs *%d ij_pairs)\n",
	 Npair_prot, N_C*(N_C+1)/2, N_ij);
  printf("Total number of classes of pairs: %d\n", *Npairs);

  /***************************** Set classes ******************************/
  int npairs=0, nsite=0;
  for(int i_L=0; i_L<N_L; i_L++){
    for(int i_U=0; i_U<N_U; i_U++){
      int ip=i_L*N_U+i_U;
      struct prot_class *prot=prot_class+ip;
      prot->i_L=i_L;
      prot->i_U=i_U;
      Initialize_prot_class(prot);
      prot->site=(*sites)+nsite; nsite+=N_C;
      struct site *site=prot->site;
      prot->pairs=(*pairs)+npairs; npairs+=Npair_prot;
      struct pair_class *pair=prot->pairs;
 
      for(int i_c=0; i_c<N_C; i_c++){
	Initialize_site(site, prot, i_c);
	site++;
	for(int j_c=i_c; j_c<N_C; j_c++){
	  for(int C_nat=0; C_nat<3; C_nat++){
	    for(int ij=0; ij<N_ij; ij++){
	      Initialize_pair(pair, prot);
	      pair->ij_index=ij;
	      pair->C_nat=C_nat;
	      pair->i1=prot->site+i_c;
	      pair->i2=prot->site+j_c;
	      //pair->i_c=i_c;
	      //pair->j_c=j_c;
	      pair++;
	    } // close ij
	  } // close j_c
	} // close i_c
      } // close C_nat
    } // close i_U
  } // close i_L


  /******************** Sample proteins ************************/
  double Npairs_all=0;
  int **C_mat=Allocate_mat2_i(L_max_act, L_max_act);
  int **C2=Allocate_mat2_i(L_max_act, L_max_act);
  /*int **k_indirect=Allocate_mat2_i(L_max, L_max);
    int **l_indirect=Allocate_mat2_i(L_max, L_max);*/
  for(int i_pdb=0; i_pdb<N_pdb; i_pdb++){
    struct protein *prot=prot_ptr+i_pdb;
    if((prot->nres<L_min)||(prot->nres>L_max))continue;
    int L=prot->nres;
    //if((prot->U_ave<U_min)||(prot->U_ave>U_max))continue;

    int p_class=Find_p_class(L, prot->U_ave, N_L, L_bin, N_U, U_bin);
    struct prot_class *prot_c=prot_class+p_class;

    Fill_C_mat(C_mat, prot->nres, prot->cont_list, prot->Nc);
    Fill_C_mat(C2, prot->nres, prot->indirect_contacts, prot->Nc_indirect);
    Npairs_all+=Sum_pairs(prot_c, Npair_prot, prot, L,
			  C_mat, C2, N_ij, ij_bin, IJ_MIN, N_C, NULL);

    /*Indirect_contacts(C2, k_indirect, l_indirect, &(prot->indirect_contacts),
      &prot->Nc_indirect, prot->nres, C_mat); */
  }

  printf("Sampling %.0f pairs in the PDB\n", Npairs_all);
  printf("%d length classes, L= < %d - > %d\n", N_L, L_short, L_long);
  printf("ij: %d %d Nij= %d\n", IJ_MIN, ij_bin[N_ij-1], N_ij);
  printf("%d contact energy classes\n", N_U);
  printf("%d protein classes, %d site and %d pair classes per prot\n",
	 *Np_class, N_C, Npair_prot);

  return(prot_class);
}

void Print_pairs(char *name_all, struct prot_class *prot_class, int Np_class,
		 int Npair_prot, int N_C, int PDB, float score_opt)
{
  // Output
  char name_file[400], tmp[80];

  int C_nat=0;
  if(prot_class->pairs->C_nat >=0)C_nat=1; // Native contacts are known

  // Score: 0=mut_inf 1=mut_inf*P_C_nat 2=d_lik 3=P_C_nat
  float score[Npair_prot]; int i_rank[Npair_prot];
  int nscore=4, iscore=3; char *name_score[nscore], rank_score[20];
  for(int i=0; i<nscore; i++)name_score[i]=malloc(20*sizeof(char));
  strcpy(name_score[0], "mut_inf");
  strcpy(name_score[1], "mut_inf*P_C_nat");
  strcpy(name_score[2], "d_lik");
  strcpy(name_score[3], "P_C_nat");
  strcpy(rank_score,name_score[iscore]);

  char head[1000];
  sprintf(head, "# pairs <= %.2f*%d samples are excluded\n", FTHR, Na2);
  strcat(head,"# For optimizing Lambda, each pair is weighted with ");
  if(I_WEIGHT==0){strcat(head, "uniform weight\n");}
  else if(I_WEIGHT==1){strcat(head, "number of pairs\n");}
  else if(I_WEIGHT==2){strcat(head, "sqrt(n)\n");}
  else if(I_WEIGHT==3){strcat(head, "mutual information\n");}
  strcat(head, "# Final scores weighted with number of pairs\n");
  strcat(head,"# Global Q matrix is ");
  if(Q_GLOB){strcat(head," used\n");}
  else{strcat(head," not used\n");}

  /*********************************************************
         Properties of prot_classes in
     <input>_NC<>_NL<>_NU<>_QGLOB<>_<misfold>_<score>_prots.dat
  ************************************************************/
  struct prot_class *prot;
  FILE *file_out; int i, j, c;
  char string[100000];
  float d_thr=8;
  int tot_predict=0, N_pred[3], N_true[3];

  // Write properties of prot_classes: Lambda, dKL, r^2, MSE
  sprintf(name_file, "%s_prots.dat", name_all);
  file_out=fopen(name_file, "w");
  fprintf(file_out, "%s", head);
  fprintf(file_out, "###1=C_nat\t2=lik\t3=<dKL>\t4=<MSE>"
	  "\t5=<r_obs_pred>\t6=<r_obs_ene>\t7=<cont_score>\t8=Npair\n");

  int N_pair_all[3], N_pred_all[3], N_true_all[3];
  for(c=0; c<3; c++){
    N_pair_all[c]=0; N_pred_all[c]=0; N_true_all[c]=0;
  }
  for(i=0; i<Np_class; i++){
    prot=prot_class+i; if(prot->Np==0)continue;
    printf("Protein class %d, L=%.1f U=%.2g Lambda=%.2f\n",
	   i, prot->L, prot->U_ave, prot->Lambda);
    struct pair_class *pair= prot->pairs;

    // Predicted contacts
    for(c=0; c<3; c++){
      prot->Npair[c]=0; N_pred[c]=0; N_true[c]=0;
    }

    if(PDB){ // PDB
      for(j=0; j<Npair_prot; j++){
	if(pair->n<=NTHR)continue;
	prot->Npair[pair->C_nat]++;
	if(pair->C_pred<0 && pair->C_pred>2){
	  printf("ERROR, c_pred= %d\n",pair->C_pred);
	}else{
	  N_pred[pair->C_pred]++;
	  if(pair->C_pred==pair->C_nat)N_true[pair->C_nat]++;
	}
	pair++;
      }
    }else{ // MSA
      int N_predict=Nc1L[(int)prot->L]*Pred_frac;
      if(N_predict>Npair_prot)N_predict=Npair_prot;
      Load_score(score, prot->pairs, Npair_prot, iscore);
      Sort(i_rank, score, Npair_prot);
      if(C_nat){
	for(j=0; j<Npair_prot; j++){
	  struct pair_class *pp=pair+i_rank[j];
	  float d=pp->min_dist;
	  if(d<=d_thr){
	    prot->Npair[1]++;
	    if(j<N_predict)N_true[1]++;
	  }else{
	    prot->Npair[0]++;
	    if(j>=N_predict)N_true[0]++;
	  }
	}
      }
      tot_predict+=N_predict;
      N_pred[1]=N_predict; N_pred[0]=Npair_prot-N_predict;
    } // end MSA

    fprintf(file_out,
	    "#L= %.0f <U>= %.3f Lambda= %.3f f_indir= %.3g Np= %d\n",
	    prot->L, prot->U_ave*TEMP, prot->Lambda, prot->f_indir, prot->Np);
    for(c=0; c<3; c++){
      fprintf(file_out, "%d\t%.5f\t%.5f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\n", c,
	      prot->lik[c], prot->dKL[c], prot->MSE[c],
	      prot->r_op[c], prot->r_oe[c], prot->cont_score[c],
	      prot->Npair[c]);
    }
    for(c=0; c<3; c++){
      N_pair_all[c]+=prot->Npair[c];
      N_pred_all[c]+=N_pred[c];
      N_true_all[c]+=N_true[c];
      sprintf(string, "# C_nat= %d pred %d", c, N_pred[c]);
      if(C_nat){
	sprintf(tmp, " true= %d positives %d", N_true[c], prot->Npair[c]);
	strcat(string, tmp);
	sprintf(tmp, " true/pred= %.3f", (float)N_true[c]/N_pred[c]);
	strcat(string, tmp);
	sprintf(tmp, " true/positive= %.3f",(float)N_true[c]/prot->Npair[c]);
	strcat(string, tmp);
      }
      printf("%s\n", string);
      fprintf(file_out, "%s\n", string);
    }
  }
  fclose(file_out);
  printf("Writing prot classes in file  %s\n", name_file);

  /*********************************************************
         Average properties over all prots in
     <input>_NC<>_NL<>_NU<>_QGLOB<>_<misfold>_<score>_ave.dat
  ************************************************************/
  // Compute and write average properties for each protein: dKL, lik
  sprintf(name_file, "%s_ave.dat", name_all);
  file_out=fopen(name_file, "w");
  fprintf(file_out, "%s", head);
  printf("Writing average properties in file  %s\n", name_file);
  //fprintf(file_out, "# %d predicted pairs with d<= %.2f based on %s\n",
  //	  tot_predict, d_thr, rank_score);
  double lik=0, dKL=0, MSE=0, norm=0, cont_score=0, norm1=0;
  int Np=0;
  for(i=0; i<Np_class; i++){
    prot=prot_class+i;
    Np+=prot->Np;
    cont_score+=prot->cont_score[1]*prot->norm[1];
    norm1+=prot->norm[1];
    for(int c=0; c<3; c++){
      MSE+=prot->MSE[c]*prot->norm[c];
      lik+=prot->lik[c]*prot->norm[c];
      dKL+=prot->dKL[c]*prot->norm[c];
      norm+=prot->norm[c];
    }
  }
  //fprintf(file_out, "#I_WEIGHT=%d I_SCORE=%d %d pairs per prot\n",
  //	  I_WEIGHT, I_SCORE, Npair_prot);
  fprintf(file_out, "#REM\tTEMP\t<lik>\t<dKL>\t<MSE>\tproteins\tnorm");
  if(C_nat){
    fprintf(file_out, "\ttrue_pred_cont\twrong_pred_cont\twrong_pred\tC_score");
  }
  fprintf(file_out, "\n");
  fprintf(file_out, "%d\t%.2f\t%.6f\t%.6f\t%.4f\t%d\t%.3g",
	  REM,TEMP,lik/norm,dKL/norm,MSE/norm,Np,norm);
  if(C_nat){
    float true_pos=N_true_all[1], false_pos=N_pred_all[1]-N_true_all[1];
    if(N_pair_all[1]){true_pos/=N_pair_all[1];}
    if(N_pred_all[1]){false_pos/=N_pred_all[1];}
    fprintf(file_out,"\t%.3f\t%.3f",true_pos,false_pos);
    int trues=0, all=0;
    for(c=0; c<3; c++){
      all+=N_pair_all[c];
      trues+=N_true_all[c];
    }
    fprintf(file_out,"\t%.3f", 1-(float)trues/all);
    fprintf(file_out,"\t%.4f", cont_score/norm1);
  }
  fprintf(file_out,"\n");
  fclose(file_out);

  /*********************************************************
         Properties of pair_classes (for each protein) in
     <input>_NC<>_NL<>_NU<>_QGLOB<>_<misfold>_<score>_pairs.dat
  ************************************************************/
  // Write properties of pair_classes: dKL, r^2, MSE, P(cont|pairs)
  FILE *file_pair[3], *file_ali;
  char header[1000]; int k;
  sprintf(header, "###"); 
  if(PDB){
    k=1; //strcat(header, "1=L\t2=[U]\t3=Lambda");
  }else{
    k=4; strcat(header,"rank\ti\tj\t"); // Alignment
  }
  sprintf(tmp, "%d=nc1\t%d=nc2\t%d=|i-j|", k, k+1, k+2); k+=3;
  strcat(header, tmp);
  sprintf(tmp, "\t%d=<Cij>\t%d=<C*Nc>-<C><Nc>\t%d=<ni*nj>-<ni><nj>",
	  k, k+1, k+2); k+=3;
  strcat(header, tmp);
  //sprintf(tmp,   "\t%d=<C*ni>-<C><ni>", ++k); //\t13=err_conv
  
  sprintf(tmp, "\t%d=E_cont\t%d=mut_inf", k, k+1); k+=2;
  strcat(header, tmp);
  sprintf(tmp, "\t%d=lik\t%d=dKL\t%d=MSE", k, k+1, k+2); k+=3; //offset
  strcat(header, tmp);
  sprintf(tmp, "\t%d=r_obs_pred\t%d=slope_op", k, k+1); k+=2;
  strcat(header, tmp);
  sprintf(tmp, "\t%d=r_obs_ene\t%d=slope_oe", k, k+1); k+=2;
  strcat(header, tmp);
  //strcat(header, \t17=logP(C_pred|MSA)
  sprintf(tmp, "\t%d=n", k); k++;
  strcat(header, tmp);
  sprintf(tmp, "\t%d=lik[C=1]-lik[C=0]\t%d=C_score\t%d=C_pred",k,k+1,k+2);
  k+=3;
  strcat(header, tmp);
  if(PDB==0 && C_nat){
    sprintf(tmp, "\t%d=C_nat\t%d=min_dist", k, k+1); k+=2;
    strcat(header, tmp);
  }
  int lines=k;
  if(PDB){
    for(i=0; i<3; i++){
      sprintf(name_file, "%s_pairs_Cnat%d.dat", name_all, i);
      file_pair[i]=fopen(name_file, "w");
      printf("Writing pair classes with Cnat= %d in file %s\n",i, name_file);
      fprintf(file_pair[i], "%s", head);
      fprintf(file_pair[i], "%s\n", header);
    }
  }else{
    // Alignment
    sprintf(name_file, "%s_pairs.dat", name_all);
    file_ali=fopen(name_file, "w"); file_out=file_ali;
    fprintf(file_out, "%s", head);
    printf("Writing pair classes in file  %s\n", name_file);
    fprintf(file_out, "%s\n", header); 
    fprintf(file_out,"###0 0 0 ");
    for(i=3; i<lines; i++)fprintf(file_out," 1");
    if(C_nat){fprintf(file_out," 0 1");}
    fprintf(file_out,"\n");
  }

  k=0;
  for(i=0; i<Np_class; i++){
    prot=prot_class+i;
    if(prot->Np==0){continue;}
    if(PDB==0){
      Load_score(score, prot->pairs, Npair_prot, iscore);
      //Rank(i_rank, score, Npair_prot);
      Sort(i_rank, score, Npair_prot);
      fprintf(file_out,"#Pairs sorted according to score %s\n", rank_score);
      fprintf(file_out,"# Prot_class %d L= %.0f [U]= %.3f "
	      "Lambda=%.3f %d proteins\n",
	      i, prot->L, prot->U_ave*TEMP, prot->Lambda, prot->Np);
    }else{
      for(j=0; j<3; j++){
	file_out=file_pair[j];
	fprintf(file_out,"# Prot_class %d L= %.0f [U]= %.3f "
		"Lambda=%.3f %d proteins\n",
		i, prot->L, prot->U_ave*TEMP, prot->Lambda, prot->Np);
      }
    }
    struct pair_class *pair;
    for(j=0; j<Npair_prot; j++){
      if(PDB){
	pair=prot->pairs+j; if(pair->n<=NTHR)continue;
	file_out=file_pair[pair->C_nat];
	//fprintf(file_out, "%d\t",k);
	//fprintf(file_out,"%.0f\t%.3f\t%.2f", //\t%.2g
	//	prot->L, prot->U_ave*TEMP, prot->Lambda);//, prot->f_indir
      }else{
	pair= prot->pairs+i_rank[j];
	if(pair->n<=NTHR)continue;
	fprintf(file_out, "%d\t%d\t%d\t",
		i_rank[j], pair->i1->i, pair->i2->i);
      }
      fprintf(file_out, "%.1f\t%.1f\t%.3g",
	      pair->i1->nc, pair->i2->nc, pair->ij);
      fprintf(file_out, "\t%.4f\t%.3f\t%.3f",
	      pair->Cont_freq, pair->Cont_freq_Nc, pair->nc_nc);
      //fprintf(file_out, "\t%.2g", pair->Cont_freq_nc);
      //fprintf(file_out, "\t%.2g", pair->error);
      fprintf(file_out,"\t%.4f", pair->E_cont_obs);
      fprintf(file_out,"\t%.4f", pair->mut_inf);

      int C; if(C_nat){C=pair->C_nat;}else{C=pair->C_pred;}
      fprintf(file_out, "\t%.4f\t%.4f\t%.3g",
	      pair->score.log_lik[C],
	      pair->score.dKL[C],
	      pair->score.MSE[C]);
      fprintf(file_out,"\t%.3f\t%.2f",
	      pair->score.r_op, pair->score.slope_op);
      fprintf(file_out,"\t%.3f\t%.2f",
	      pair->score.r_oe, pair->score.slope_oe);
      //pair->score.offset, 
      fprintf(file_out,"\t%8.0f", pair->n);
      float lP=log(pair->score.P_C_nat);
      if((lP<-1000) ||isnan(lP))lP=-1000;
      //fprintf(file_out,"\t%.2g", lP);
      fprintf(file_out,"\t%.4f", pair->score.log_lik[1]-pair->score.log_lik[0]);
      fprintf(file_out,"\t%.4f", pair->cont_score);
      fprintf(file_out,"\t%d", pair->C_pred);
      if((PDB==0)&&(C_nat))
	fprintf(file_out,"\t%d\t%.1f", pair->C_nat, pair->min_dist);
      fprintf(file_out,"\n");
      k++;
    }
  }
  if(PDB){
    for(i=0; i<3; i++)fclose(file_pair[i]);
  }else{
    fclose(file_out);
  }

  if(PRINT_Q){
    // Print unspecific Q
    printf("Writing global Q_obs in file ");
    if(ALL_Q){
      sprintf(name_file, "%s_Q_global.dat", name_all);
      printf("%s\n", name_file);
      FILE *file_out=fopen(name_file, "w");
      fprintf(file_out,"#log_Q_global\n");
      Print_log_Q(NULL, prot_class->log_Q_glob, file_out, 0);
      fclose(file_out);
    }else{
      printf("%s_Q_global_Protk.dat 0<=k<%d\n", name_all, Np_class);
      for(i=0; i<Np_class; i++){
	prot=prot_class+i;
	if(prot->Np==0){continue;}
	sprintf(name_file, "%s_Q_global_Prot%d.dat", name_all, i);
	FILE *file_out=fopen(name_file, "w");
	fprintf(file_out,"#L= %.0f [U]= %.3f Lambda=%.3f %d proteins\n#Qobs\n",
		prot->L, prot->U_ave*TEMP, prot->Lambda, prot->Np);
	Print_log_Q(NULL, prot->log_Q_glob, file_out, 0);
	fclose(file_out);
      }
    }
  
    if(PDB){
      // Write Q_obs and Q_pred for each class
      printf("Writing Q_obs and Q_pred in file  %s_Q_<>.dat", name_all);
      printf("   0<= k < %d\n", Np_class*Npair_prot);  k=0;
      for(i=0; i<Np_class; i++){
	prot=prot_class+i;
	if(prot->Np==0){continue;}
	for(int C_nat=0; C_nat<3; C_nat++){
	  for(int i_c=0; i_c<N_C; i_c++){
	    for(int j_c=i_c; j_c<N_C; j_c++){
	      sprintf(name_file, "%s_Q_Prot%d_nc%d%d_Cnat%d.dat", name_all,
		      i, i_c, j_c, C_nat);
	      FILE *file_out=fopen(name_file, "w");
	      fprintf(file_out,"#L= %.0f [U]= %.3f Lambda=%.3f %d proteins\n",
		      prot->L, prot->U_ave*TEMP, prot->Lambda, prot->Np);
	      fprintf(file_out,"#C_nat= %d\n", C_nat);
	      for(j=0; j<Npair_prot; j++){
		struct pair_class *pair= prot->pairs+j;
		if(pair->n<=NTHR || pair->C_nat!=C_nat ||
		   pair->i1->i_c!= i_c || pair->i2->i_c!= j_c){
		  continue;
		}
		fprintf(file_out,
			"# ij=%.0f n= %.0f C_pred= %d r_op=%.2f MSE=%.2f\n",
			pair->ij, pair->n, pair->C_pred,
			pair->score.r_op,pair->score.MSE[pair->C_pred]);
	      }
	      fprintf(file_out, "#AA1AA2\t2=Econt "); k=3;
	      for(j=0; j<Npair_prot; j++){
		struct pair_class *pair= prot->pairs+j;
		if(pair->n<=NTHR || pair->C_nat!=C_nat ||
		   pair->i1->i_c!= i_c || pair->i2->i_c!= j_c){
		  continue;
		}
		fprintf(file_out,
			"\t%d=E_not_ave_ij%.0f\t"
			"%d=Qobs_ij%.0f\t%d=Qpred_ij%.0f",
		      k, pair->ij, k+1, pair->ij, k+2, pair->ij); k+=3;
	      }
	      fprintf(file_out, "\n");
	      for(int ab=0; ab<Na2; ab++){
		int a=ab/Naa, b=ab-(a*Naa);
		fprintf(file_out, "%c%c\t%.3f\t",
			AMIN_CODE[a],AMIN_CODE[b], E_cont_T[a][b]);
		for(j=0; j<Npair_prot; j++){
		  struct pair_class *pair= prot->pairs+j;
		  if(pair->n<=NTHR || pair->C_nat!=C_nat ||
		     pair->i1->i_c!= i_c || pair->i2->i_c!= j_c){
		    continue;
		  }
		  float E_cont_not[Na2];
		  Zero_mean(E_cont_not, E_cont_1,
			    pair->P1i1_obs, pair->P1i2_obs, Naa);
		  fprintf(file_out, "%.3f\t%.3f\t%.3f\t",
			  E_cont_not[ab],
			  pair->log_Q_obs[ab],
			  pair->log_Q_pred[pair->C_pred][ab]);
		} // end pairs
		fprintf(file_out,"\n");
	      } // end ab
	      fclose(file_out);  

	    } // end j_c
	  } // end i_c
	} // end C_nat
      } // end prot
    } // end if(PDB)
  } // end if(PRINT_Q)

  if(C_nat){
    // Contact overlap
    sprintf(name_file, "%s_contact_overlap.dat", name_all);
    file_out=fopen(name_file, "w");
    printf("Writing contact overlap and MCC in file  %s\n", name_file);
    int n=Npair_prot;
    for(i=0; i<Np_class; i++){
      prot=prot_class+i;
      if(prot->Np==0){continue;}
      fprintf(file_out, "# L=%.0f <U>=%.3f Lambda=%.3f Np=%d\n",
	      prot->L, prot->U_ave*TEMP, prot->Lambda, prot->Np);
      fprintf(file_out, "# Pairs ranked according to score %s\n",
	      name_score[iscore]);
      Load_score(score, prot->pairs, Npair_prot, iscore);
      float max=score[0], min=score[0];
      for(j=1; j<Npair_prot; j++){
	if(score[j]>max)max=score[j];
	if(score[j]<min)min=score[j];
      }
      min*=0.5;
      float delta=(max-min)/50;

      for(float dthr=4.5; dthr< 9; dthr+=1.0){
	fprintf(file_out, "# Distance threshold: %.3f\n", dthr);
	fprintf(file_out, "# thr overlap MCC true_pos pred_cont\n");
	int n1=0; struct pair_class *p=prot->pairs;
	for(j=0; j<Npair_prot; j++){if(p->min_dist<dthr)n1++; p++;}
	fprintf(file_out, "# %d contacts %d pairs\n", n1, n);
	float max_MCC=0, xm=-1; int m=0;
	for(float x=min; x<max; x+=delta){
	  int q1=0, q0=0, m1=0; p=prot->pairs;
	  for(j=0; j<Npair_prot; j++){
	    if(score[j]>=x){m1++; if(p->min_dist<dthr)q1++;}
	    else{if(p->min_dist>=dthr)q0++;}
	    p++;
	  }
	  if((m1==m)||(m1==n)){continue;}
	  m=m1;
	  float norm=(float)m1*n1*(n-m1)*(n-n1);
	  float MCC=((q1*q0)-((m1-q1)*(n-m1-q0)))/sqrt(norm);
	  fprintf(file_out, "%.4f %.3f %.4f %d %d\n",
		  x, q1/sqrt(n1*m1), MCC, q1, m1);
	  if(MCC > max_MCC){max_MCC=MCC; xm=x;}
	}
	fprintf(file_out, "&\n");
	if(PDB==0)
	  printf("Max value of MCC= %.3f dthr= %.2f score %.3g\n",
	       max_MCC, dthr, xm);
      }
    }
    fclose(file_out);
  }

  if(C_nat && (PDB==0)){
    // Distributions
    float d0=4.5, dmax=8.5, dstep=1.00;
    int Nbin=(dmax-d0)/dstep+1;
    int numbin[Nbin+1];
    double lij[Nbin+1], lij2[Nbin+1];

    sprintf(name_file, "%s_predicted_contacts.dat", name_all);
    file_out=fopen(name_file, "w");
    printf("Writing predicted contacts in file %s\n", name_file);
    fprintf(file_out, "#Sequence identity: min= %.3f max= %.3f\n",
	    SI_min, SI_max);
    for(i=0; i<Np_class; i++){
      prot=prot_class+i;
      if(prot->Np==0){continue;}
      struct pair_class *pair=prot->pairs;
      int N_pred=Nc1L[(int)prot->L]*Pred_frac;
      if(N_pred>Npair_prot)N_pred=Npair_prot;
      fprintf(file_out, "#Prot%d L=%.0f <U>=%.3f Lambda=%.3f %d proteins",
	      i, prot->L, prot->U_ave*TEMP, prot->Lambda, prot->Np);
      fprintf(file_out, " score= %.4g\n", score_opt);
      int n1=0; struct pair_class *p=pair;
      for(j=0; j<Npair_prot; j++){if(p->min_dist<d0)n1++; p++;}
      fprintf(file_out, "#%d predictions %d contacts %.0f sequences\n",
	      N_pred, n1, pair->n);
      fprintf(file_out, "# Only specific likelihood is scored\n"); 
      for(int is=0; is<=3; is++){
	Load_score(score, prot->pairs, Npair_prot, is);
	Sort(i_rank, score, Npair_prot);
	fprintf(file_out, "# score %d: %s\n", is, name_score[is]);
	fprintf(file_out, "#d num(D<d) frac ave(|i-j|) sigma(|i-j|)\n");
	for(j=0; j<=Nbin; j++){
	  numbin[j]=0; lij[j]=0; lij2[j]=0;
	}
	double dave=0; float d;
	for(j=0; j<N_pred; j++){
	  d=(pair+i_rank[j])->min_dist;
	  k= (d-d0)/dstep;
	  if(k<0){k=0;}
	  else if(k>=Nbin){k=Nbin; dave+=d;}
	  numbin[k]++;
	  int l=(pair+i_rank[j])->ij;
	  lij[k]+=l; lij2[k]+=l*l; 
	}
	d=d0;
	for(j=0; j<=Nbin; j++){
	  if(numbin[j]){
	    lij[j]/=numbin[j]; lij2[j]-=lij[j]*lij[j];
	    if(numbin[j]>1)lij2[j]=sqrt(lij2[j]/(numbin[j]-1));
	  }
	  fprintf(file_out, "%.2f\t%d\t%.3f\t%.1f\t%.1f\n",
		  d, numbin[j], (float)numbin[j]/N_pred, lij[j], lij2[j]);
	  d+=dstep;
	  if(d>dmax){d=dave; if(numbin[Nbin])d/=numbin[Nbin];}
	}
	fprintf(file_out, "&\n");
      } // end score
    }
  }

}

void Print_average_Q(char *name_all, struct prot_class *prot_class,
		     int Np_class, int Npair_prot)
{
  // Write Q_obs and Q_pred for each class
  printf("Writing Q_obs and Q_pred in file  %s_Q_k.dat", name_all);
  printf("   0<= k < %d\n", Np_class*N_ij*2);
  int k=0, c;
  float Q_pred_all[Na2], Q_obs[Na2];
  for(int ij=0; ij<N_ij; ij++){
    for(int contact=0; contact<=1; contact++){
      double norm=0, l=0; int n=0;
      for(c=0; c<Na2; c++){Q_pred_all[c]=0; Q_obs[c]=0;}
      for(int ip=0; ip<Np_class; ip++){
	struct prot_class *prot=prot_class+ip;
	if(prot->Np==0)continue;
	for(int j=0; j<Npair_prot; j++){
	  struct pair_class *pair= prot->pairs+j;
	  if(pair->C_nat!=contact)continue;
	  if(pair->n==0)continue;
	  float w=pair->n; // pair->w;
	  int lij=Find_bin_i(pair->ij, N_ij, ij_bin);
	  if(lij!=ij)continue;
	  for(c=0; c<Na2; c++){
	    Q_pred_all[c]+=w*exp(pair->log_Q_pred[pair->C_pred][c]);
	    Q_obs[c]+=w*pair->Q_obs[c];
	  }
	  norm+=w; l+=w*pair->ij; n++;
	}
      } // end sum
      l/=norm;
      for(c=0; c<Na2; c++){
	Q_pred_all[c]/=norm;
	Q_obs[c]/=norm;
      }
      char name_file[80];
      sprintf(name_file, "%s_Q_%d.dat", name_all, k);
      FILE *file_out=fopen(name_file, "w");
      fprintf(file_out, "# <ij>= %.1f C_nat=%d %d pairs\n", l, contact, n);
      fprintf(file_out, "#Qpred Qobs\n");
      int i_rank[Na2];
      Sort(i_rank, Q_obs, Na2);
      for(int i=0; i<Na2; i++){
	int ab=i_rank[i], a=ab/Naa, b=ab-(a*Naa);
	fprintf(file_out, "%.3f\t%.3f\t%c%c\n",
		Q_pred_all[ab],Q_obs[ab],
		AMIN_CODE[a],AMIN_CODE[b]);
      }
      fclose(file_out);
      k++;
    } // end contact
  } // end ij
}


void Print_log_Q(float *log_Q_pred, float *log_Q_obs, FILE *file_out, int symm)
{
  int i_rank[Na2];
  //Rank(i_rank, log_Q_obs, Na2);
  Sort(i_rank, log_Q_obs, Na2);
  for(int i=0; i<Na2; i++){
    int ab=i_rank[i], a=ab/Naa, b=ab-(a*Naa);
    if((symm==0)&&(b>a))continue;
    if(log_Q_pred)fprintf(file_out, "%.3g\t", log_Q_pred[ab]);
    fprintf(file_out, "%.4g\t%c%c\n",log_Q_obs[ab],AMIN_CODE[a],AMIN_CODE[b]);
  }
}

void Load_score(float *score, struct pair_class *pair, int Npair_prot, int type)
{
  struct pair_class *p=pair;
  for(int j=0; j<Npair_prot; j++){
    float s;
    if(type==0){ // mut_inf
      s=p->mut_inf;
    }else if(type==1){ // mut_inf P_C_nat
      s=p->mut_inf*p->score.P_C_nat;
    }else if(type==2){ // d_lik=log_lik[1]-log_lik[0]
      s=p->score.log_lik[1]-p->score.log_lik[0];
    }else if(type==3){ // P_C_nat
      s=p->score.P_C_nat;
      /*}else if(type==6){ // Consider r
      s=p->score.P_C_nat*(p->mut_inf);
      s*=p->prot_class->Cont_freq[(int)p->ij];
      s*=p->score.r; */
    }

    /*if(p->ij<ij_max){
      score[j]*=p->prot_class->Cont_freq[(int)p->ij];
    }else{
      score[j]*=p->prot_class->Cont_freq[ij_max];
      }*/
    score[j]=s; p++;
  }
}

float Optimize_Lambda(float *La_opt, struct prot_class *prot,
		      int Npair_prot, int i_prot)
{
  /**********************  Inform  **********************/
  printf("i_prot=%d %d pairs, ",i_prot, prot->npairs);
  if(I_SCORE==0){printf("optimizing likelihood\n");}
  else if(I_SCORE==1){printf("optimizing dKL(Q_pred, Q_obs)\n");}
  else if(I_SCORE==2){printf("optimizing MSE(Q_pred, Q_obs)\n");}
  else{printf("ERROR, wrong score for optimization\n"); exit(8);}
  //printf("Maximum allowed value of Lambda: %.2f\n", LAMBDA_MAX);


  /**********************  Optimize Lambda  **********************/
  // Parameters
  int IT_MAX=20, iter, n_fail=0, n_fail_max=1; //float EPS=0.0005;
  float La_STEP=0.2, La_inf=0.01, La_sup=LAMBDA_MAX, LAMBDA=1.4;

  float Lambda=LAMBDA, y=Score(Lambda, prot, Npair_prot, I_SCORE, 0);
  float L0=Lambda, y0=y, y_opt=y; *La_opt=L0;
  printf("it=0 Lambda=%.3f Score=%.6g\n", Lambda, y);  

  Lambda+=La_STEP;
  y=Score(Lambda, prot, Npair_prot, I_SCORE, 0);
  float L1=Lambda, y1=y; if(y > y_opt){*La_opt=Lambda; y_opt=y;}
  printf("it=0 Lambda=%.3f Score=%.6g\n", Lambda, y);  

  Lambda+=La_STEP;
  y=Score(Lambda, prot, Npair_prot, I_SCORE, 0);
  float L2=Lambda, y2=y; if(y > y_opt){*La_opt=Lambda; y_opt=y;}
  printf("it=0 Lambda=%.3f Score=%.6g\n", Lambda, y);  

  for(iter=0; iter<IT_MAX; iter++){
    Lambda=Find_max_quad(L0, L1, L2, y0, y1, y2, La_inf, La_sup);
    if(Lambda==L0 || Lambda==L1 || Lambda==L2){
      printf("Lambda=%.3f repeats, leaving the loop\n", Lambda); break;
    }
    if(isnan(Lambda)||Lambda<La_inf ||Lambda>La_sup)break;
    printf("it=%d computing Lambda=%.3f\n", iter+1,Lambda);  
    y=Score(Lambda, prot, Npair_prot, I_SCORE, 0);
    printf("it=%d Lambda=%.3f Score=%.6g\n", iter+1,Lambda,y);  
    if(y > y_opt){*La_opt=Lambda; y_opt=y; n_fail=0;}
    else{n_fail++; if(n_fail >= n_fail_max){break;}}
    //if((Lambda==La_inf)||(Lambda==La_sup))break;

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
    printf("WARNING, Lambda optimization did not converge (prot %d)\n",i_prot);


  printf("### prot %d Optimal Lambda: %.3f Score= %.6g %d proteins",
	 i_prot, *La_opt, y_opt, prot->Np);
  //y_opt=Score(*La_opt, prot, Npair_prot, I_SCORE, 1);
  //printf(" confirm score= %.6g", y_opt);
  printf("\n");

  return(y_opt);
}

