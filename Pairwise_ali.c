#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "protein4.h"
#include "pairwise.h"
#include "ali_pdb.h"
#include "alignments.h"
#include "NeedlemanWunsch.h"
#include "Pairwise_aux.h"
#include "codes.h"
#include "allocate.h"
#include "REM_new.h"
#include "externals.h"

extern void Get_name(char *name, char *file_name);

int Align_fast(int *ali, char *seq_pdb, int L_pdb, char *msa, int L_ali,
	       int LMIN);
int Align_dynpro(int *ali, char *seq_pdb, int L_pdb, char *msa, int L_ali,
		 int LMIN);
char *Remove_gaps(int *L_seq, int **ali_nogap, char *msa, int L_ali);
void Restore_gaps(int *ali, int L, int *ali_nogap);
int Mut_min(char *Seq_i, char *Seq_j, int L_ali, int mutmin);
float Seq_id(char *Seq_i, char *Seq_j, int L_ali);
void Copy_seq(char *Seq_k, char *Seq_i, int L_ali);
int Seq_length(char *Seq, int L_ali);

int **Align_pdb(int *target, int **aligned, float ***seq_id,
		char **MSA, int n_seq, int L_ali, char *file_ali,
		struct protein *template, int N_template)
 {
   int Compute_SI=0; float *sid;
   if(N_template==0)return(NULL);
   int **ali_seq=malloc(N_template*sizeof(int *));
   *aligned=malloc(N_template*sizeof(int));
   if(Compute_SI)*seq_id=malloc(N_template*sizeof(float *));
   int id_max=0, s_max=-1;
   for(int i_pdb=0; i_pdb<N_template; i_pdb++){
      struct protein *prot=template+i_pdb;
      int L_pdb=prot->nres, i_seq, i;
      printf("Aligning %s, %d residues\n", prot->name, prot->nres);
      int *ali=malloc(L_pdb*sizeof(int));
      ali_seq[i_pdb]=ali; (*aligned)[i_pdb]=-1;
      int L=L_pdb; if(L_ali<L)L=L_ali;
      int id_min=0.9*L, id=0; // Minimum sequence identity
      /*for(i_seq=0; i_seq<n_seq; i_seq++){ 
	id=Align_fast(ali, prot->amm, L_pdb, MSA[i_seq], L_ali, id_min);
	if(id>=id_max){id_max=id; s_max=i_seq;}
	if((id>=id_min)&&((L-id)<3))goto end;
	}*/
      for(i_seq=0; i_seq<n_seq; i_seq++){ 
	id=Align_dynpro(ali, prot->amm, L_pdb, MSA[i_seq], L_ali, id_min);
	if(id>=id_max){id_max=id; s_max=i_seq;}
	if(id>=id_min){break;} break; // Align only with the first seq.
      }
      if(Compute_SI){
	(*seq_id)[i_pdb]=malloc(n_seq*sizeof(float));
	sid=(*seq_id)[i_pdb];
      }
      if(id<id_min){
	printf("WARNING, PDB %s not found in MSA, id=%d required: %d seq:\n",
	       prot->name, id, id_min);
	for(i=0; i<L_pdb; i++){printf("%c", prot->amm[i]);} printf("\n");
	//for(i_seq=0; i_seq<n_seq; i_seq++)sid[i_seq]=0;
      }else{
	if(s_max>=0)*target=s_max; 
	(*aligned)[i_pdb]=s_max;
	if(Compute_SI){
	  for(i_seq=0; i_seq<n_seq; i_seq++){
	    int id=0;
	    for(int j=0; j<L_pdb; j++){
	      if((ali[j]>=0)&&(MSA[i_seq][ali[j]]==prot->amm[j]))id++;
	    }
	    sid[i_seq]=(float)id/(float)L_pdb;
	  }
	}
      }
   }
   return(ali_seq);
 }

int Align_fast(int *ali, char *seq_pdb, int L_pdb, char *msa, int L_ali,
	       int LMIN)
{
  int id=0, SHIMAX=10, GMAX=10;
  int L_seq, *ali_nogap;
  char *seq=Remove_gaps(&L_seq, &ali_nogap, msa, L_ali);

  int i, j, shift; // Only gaps in sequence
  for(shift=0; shift <= SHIMAX; shift ++){
    int gap=0; id=0; j=0;
    for(i=0; i<L_pdb; i++)ali[i]=-1;
    for(i=shift; i<L_pdb; i++){   
      while(seq[j]!=seq_pdb[i]){
	j++; gap++; if((j>=L_seq)||(gap>GMAX))break;
      }
      ali[i]=j; id++; j++;
      if(j>=L_seq)break;
    }
    if(id >LMIN)goto end;
  }
  // Only gaps in structure
  for(shift=0; shift <= SHIMAX ; shift ++){
    for(i=0; i<L_pdb; i++)ali[i]=-1;
    int gap=0; id=0; i=0;
    for(j=shift; j<L_ali; j++){
      while(seq[j]!=seq_pdb[i]){
	i++; gap++;
	if((i>=L_pdb)||(gap>GMAX))break;
      }
      ali[i]=j; id++; i++;
      if(i>=L_pdb)break;
    }
    if(id > LMIN)goto end;
  }
  // Mismatch, no gap
  int gap=0; id=0; j=0;
  for(i=0; i<L_pdb; i++)ali[i]=-1;
  for(i=0; i<L_pdb; i++){   
    while(seq[j]!=seq_pdb[i]){
      ali[i]=j; i++; j++; gap++;
      if((i>=L_pdb)||(j>=L_seq)||(gap>GMAX))break;
    }
    ali[i]=j; id++; j++; if(j>=L_seq)break;
  }
  if(id >LMIN)goto end;

 end:
  if(id >LMIN)Restore_gaps(ali, L_pdb, ali_nogap);
  free(ali_nogap);
  free(seq);
  return(id);
}


int Align_dynpro(int *ali, char *seq_pdb, int L_pdb, char *msa, int L_ali,
		 int LMIN)
{
  int VBS=0;  // Verbose:
  int IDE=1;  // Use identity to score alignment
  int GAP=3;  // Gap opening penalty, reduced to favour better alis
  //int GAP=7;  // Gap opening penalty

  int L_seq, *ali_nogap, id=0;
  char *seq=Remove_gaps(&L_seq, &ali_nogap, msa, L_ali);

  char *ali1=malloc((L_pdb+L_seq)*sizeof(char));
  char *ali2=malloc((L_pdb+L_seq)*sizeof(char));
  int nali, al=
    alignNW(seq_pdb, L_pdb, seq, L_seq, VBS, IDE, GAP, ali1, ali2, &nali);
  if(al==0){
    printf("WARNING, NW alignment failed\n"); goto empty;
  }
  int i1=0, i2=0;
  for(int i=0; i<nali; i++){
    if((ali1[i]!='-')&&(ali2[i]!='-')){
      if(ali1[i]==ali2[i])id++;
    }
    if(ali2[i]=='-'){
      ali[i1]=-1;
    }else{
      ali[i1]=ali_nogap[i2]; i2++;
    }
    if(ali1[i]!='-')i1++;
  }

 empty:
  free(ali1); free(ali2);
  free(seq); free(ali_nogap);
  return(id);
}

char *Remove_gaps(int *L_seq, int **ali_nogap, char *msa, int L_ali)
{
  char *seq=malloc(L_ali*sizeof(char));
  *ali_nogap=malloc(L_ali*sizeof(int));
  int L=0;
  for(int j=0; j<L_ali; j++){
    if(msa[j]!='-'){(*ali_nogap)[L]=j; seq[L]=msa[j]; L++;}
  }
  *L_seq=L;
  return(seq);
}

void Restore_gaps(int *ali, int L, int *ali_nogap){
  for(int i=0; i<L; i++){
    if(ali[i]>=0)ali[i]=ali_nogap[ali[i]];
  }
}

void Convert_sequences(int **MSA_num, char **MSA, int n_seq, int L_ali){
  for(int i=0; i<n_seq; i++){
    int *num=MSA_num[i];
    char *aa=MSA[i];
    for(int j=0; j<L_ali; j++){
      if(*aa=='-'){*num=-1;}
      else{*num=Code_AA(*aa);}
      aa++; num++;
    }
  }
}

int Remove_gaps_target(char **MSA, int i_target, int n_seq, int L_ali,
		       int **ali, struct protein *prots, int N_template,
		       float thr)
{
  int nthr=n_seq*thr, mtar=0, mgap=0;
  int L=0, i, k; char *target=MSA[i_target];
  for(int j=0; j<L_ali; j++){
    int gap=0; // Count number of gaps at position j
    for(i=0; i<n_seq; i++)if((MSA[i][j]=='-')||(MSA[i][j]=='.'))gap++;
    if((target[j]!='-')&&(target[j]!='.')&&(gap<nthr)){
      if(j!=L){
	for(i=0; i<n_seq; i++)MSA[i][L]=MSA[i][j];
	for(i=0; i<N_template; i++){
	  int L_pdb=prots[i].nres;
	  for(k=0; k<L_pdb; k++)if(ali[i][k]==j){ali[i][k]=L; break;}
	}
      }
      L++;
    }else{
      if(gap>=nthr){mgap++;}else{mtar++;}
    }
  }
  printf("Eliminating %d columns that are not present in target\n", mtar);
  printf("Eliminating %d columns with > %.0f percent gaps\n", mgap,thr*100);
  return(L);
}

void Set_C_nat_ali(struct pair_class *pairs, int Nij_class,
		   struct protein *template, int N_template,
		   int **ali_seq)
{
  if((ali_seq==NULL)||(template==NULL))return;
  for(int ip=0; ip<N_template; ip++){
    struct protein *prot=template+ip;
    struct contact *cont=prot->cont_list;
    int *ali=ali_seq[ip];
    for(int i=0; i<prot->Nc; i++){
      if((ali[cont->res1]>=0)&&(ali[cont->res2]>=0)&&
	 (label_ij[ali[cont->res1]][ali[cont->res2]]>=0)){
	(pairs+label_ij[ali[cont->res1]][ali[cont->res2]])->C_nat=1;
      }
      cont++;
    }
    cont=prot->indirect_contacts;
    for(int i=0; i<prot->Nc; i++){
      if((ali[cont->res1]>=0)&&(ali[cont->res2]>=0)&&
	 (label_ij[ali[cont->res1]][ali[cont->res2]]>=0)){
	struct pair_class *pair=pairs+
	  label_ij[ali[cont->res1]][ali[cont->res2]];
	if(pair->C_nat<=0)pair->C_nat=2;
      }
      cont++;
    }
    for(int i=0; i<prot->nres; i++){
      if(ali[i]<0){continue;}
      int *lab_i=label_ij[ali[i]];
      for(int j=i+2; j<prot->nres; j++){
	if((ali[j]>=0)&&(lab_i[ali[j]]>=0)){
	  struct pair_class *pair=pairs+lab_i[ali[j]];
	  if((ip==0)||(prot->min_dist[i][j]<pair->min_dist)){
	    pair->min_dist=prot->min_dist[i][j];
	  }
	}
      }
    }
  }
  struct pair_class *pair=pairs;
  for(int i=0; i<Nij_class; i++){
    if(pair->C_nat<0)pair->C_nat=0;
    pair++;
  }
}

void Copy_MSA(char **MSA, char **MSA_0, int n_seq, int L_ali){
  for(int i=0; i<n_seq; i++){
    char *M=MSA[i], *M0=MSA_0[i];
    for(int j=0; j<L_ali; j++){*M=*M0; M++; M0++;}
  }
}

int Select_seqs(char **MSA, int *i_target,
		int n_seq, int L_ali, float *SI_target, float SI_min)
{
  // Eliminate sequences with SI<SI_min with respect to the target
  if(*i_target<0)return(n_seq);
  printf("Eliminating sequences with SI(s,target)< %.3f\n", SI_min);
  int k=0;
  for(int i=0; i<n_seq; i++){
    if((i==*i_target)||(SI_target[i]>=SI_min)){
      if(i!=k){
	Copy_seq(MSA[k], MSA[i], L_ali);
	SI_target[k]=SI_target[i];
	if(i==*i_target)*i_target=k;
      }
      k++;
    }
  }
  printf("Selecting %d seqs out of %d\n", k, n_seq);
  return(k);
}

float *Compute_SI_target(char **MSA, int n_seq, int L_ali, int i_target){
  float *SI_target=malloc(n_seq*sizeof(float));
  char *Seq_t=MSA[i_target];
  for(int i=0; i<n_seq; i++)SI_target[i]=Seq_id(MSA[i], Seq_t, L_ali);
  return(SI_target);
}

void Cluster_seqs(int *cluster, char **MSA, int ns,
		  int L_ali, float SI_max)
{
  printf("Clustering seqs with SI > %.3f\n", SI_max);
  for(int i=0; i<ns; i++)cluster[i]=i;
  if(SI_max > 1)return;

  float m_min=1.-SI_max;
  int L_seq[ns], i, k;
  for(i=0; i<ns; i++)L_seq[i]=Seq_length(MSA[i], L_ali);
  int nc[ns]; for(i=0; i<ns; i++)nc[i]=1;
  for(i=0; i<ns; i++){
    char *Seq_i=MSA[i]; int mutmin;
    for(int j=0; j<i; j++){
      if(cluster[j]==cluster[i])continue;
      if(L_seq[i]<L_seq[j]){mutmin=L_seq[i];}else{mutmin=L_seq[j];}
      mutmin*=m_min;
      if(Mut_min(Seq_i, MSA[j], L_ali, mutmin)==0){
	// Join cluster j with cluster i
	int c=cluster[i], c1=cluster[j];
	if(c1>c){int tmp=c1; c1=c; c=tmp;}
	if(nc[c]>1){
	  for(k=0; k<ns; k++)if(cluster[k]==c)cluster[k]=c1;
	}else{
	  cluster[i]=c1;
	  cluster[j]=c1;
	}
	nc[c1]+=nc[c];
	nc[c]=0;
      }
    }
  }

  k=0; int n=0;
  for(int c=0; c<ns; c++){
    if(nc[c]>0){
      int m=0; for(i=0; i<ns; i++)if(cluster[i]==c)m++;
      if(m!=nc[c]){
	printf("ERROR, %d elements instead of %d in cluster %d\n",
	       m, nc[c], k); exit(8);
      }
      k++; n+=nc[c];
    }
  }
  printf("Grouping %d seqs in %d clusters\n", ns, k);
  if(n!=ns){
    printf("ERROR, wrong number of elements, %d instead of %d\n",n,ns);
    exit(8);
  }
  return;
}

void Weight_seqs(float *w_seq, int *cluster, int ns, float *SI_target, int WSI)
{
  float SI_0=0.07; // Sequence identity of unrelated proteins
  int nele[ns], i, c, c_max=0; double SI[ns];
  for(i=0; i<ns; i++){nele[i]=0; SI[i]=0;}
  for(i=0; i<ns; i++){
    c=cluster[i];
    if(c>=ns){
      printf("ERROR, too large cluster index %d > %d\n",c,ns-1);
      exit(8);
    }
    nele[c]++;
    SI[c]+=SI_target[i]; 
    if(c>c_max)c_max=c;
  }

  int k=0, n=0;
  for(int c=0; c<=c_max; c++){
    if(nele[c]==0)continue;
    n+=nele[c]; float w=1./nele[c]; k++;
    // Old
    /*if(WSI){
      SI[c]/=nele[c]; float ww=(SI[c]-SI_0)/(1-SI_0); if(ww<0)ww=0; w*=ww;
    }
    for(i=0; i<ns; i++)if(cluster[i]==c)w_seq[i]=w; */
    if(WSI==0){
      for(i=0; i<ns; i++)if(cluster[i]==c)w_seq[i]=w;
    }else{
      for(i=0; i<ns; i++)if(cluster[i]==c){
	  w_seq[i]=w*(SI_target[i]-SI_0); if(w_seq[i]<0)w_seq[i]=0;
	}
    }
  }
  if(n!=ns){
    printf("ERROR, %d elements instead of %d in %d clusters\n",n, ns, k);
    exit(8);
  }
}

void Count_Pairwise(struct pair_class *pairs, int Npairs,
		    struct prot_class *prot_class, int **MSA_aa,
		    int n_seq, float *w_seq, int L_tar,
		    int ij_min, int Naa) //, int PURGE_REDUNDANCE
{
  int i, j;
  struct pair_class *pair=pairs;
  for(i=0; i<L_tar; i++){
    for(j=i+ij_min; j<L_tar; j++){
      pair->E_cont_obs=0;  float norm=0;
      for(int a=0; a<Na2; a++)pair->N2_obs[a]=0;
      for(int i_seq=0; i_seq<n_seq; i_seq++){
	float w=w_seq[i_seq];
	int aa_i=MSA_aa[i_seq][i], aa_j=MSA_aa[i_seq][j];
	if((aa_i>=20)||(aa_j>=20))continue;
	// gap
	if(aa_i<0){if(Naa==21){aa_i=20;}else{continue;}}
	if(aa_j<0){if(Naa==21){aa_j=20;}else{continue;}}
	int c=aa_i*Naa+aa_j;
	if((c<0)||(c>=Na2)){
	  printf("ERROR, wrong ab= %d\n", c); exit(8);
	}
	pair->N2_obs[c]+=w;
	norm+=w;
	pair->E_cont_obs+=E_cont_T[aa_i][aa_j]*w;
      }
      pair->E_cont_obs/=norm;
      pair++;
    }
  }
}

struct pair_class *
Set_pairs_ali(int *Np_class, int *Npairs, struct prot_class **prot_class,
	      int ij_min, int ij_max, int L_tar, int L, struct protein *prot)
{
  printf("Setting pairs of residues for MSA. L= %d\n", L_tar);
  printf("ij_min= %d L_max= %d Naa= %d\n", ij_min, L, Naa);

  /**************************  Allocate  ****************************/
  //Na2=Naa*Na
  
  *Np_class=1;
  *prot_class=malloc(1*sizeof(struct prot_class));
  Initialize_prot_class(*prot_class);

  struct site *sites=malloc(L_tar*sizeof(struct site));
  (*prot_class)->site=sites;

  *Npairs=0;
  label_ij=malloc(L_tar*sizeof(int *));
  for(int i=0; i<L_tar; i++){
    label_ij[i]=malloc(L_tar*sizeof(int));
    for(int j=0; j<i+ij_min; j++){
      label_ij[i][j]=-1;
      label_ij[j][i]=-1;
    }
    for(int j=i+ij_min; j<L_tar; j++){
      label_ij[i][j]=(*Npairs);
      label_ij[j][i]=(*Npairs);
      (*Npairs)++;
    }
  }
  struct pair_class *pair_class=malloc((*Npairs)*sizeof(struct pair_class));
  struct pair_class *pair=pair_class;
  for(int i=0; i<(*Npairs); i++){
    Initialize_pair(pair, (*prot_class));
    pair++;
  }


  Sum_pairs(*prot_class, *Npairs, prot, L_tar,
	    NULL, NULL, 0, NULL, ij_min, 0, label_ij);

  printf("%d classes of residue pairs\n", *Npairs);
  return(pair_class);
}


void Find_indirect(struct pair_class *pairs, int Nij_class)
{
  int a, b, i, j;
  struct pair_class *pair=pairs;
  for(a=0; a<Nij_class; a++){
    pair->P_ind=0; pair->ik=-1; pair->jk=-1; pair++;
  }
  for(a=0; a<Nij_class; a++){
    struct pair_class *pair1=pairs+a;
    if(pair1->mut_inf<0)continue; // Comment if selection on P_C_nat
    for(b=a+1; b<Nij_class; b++){
      struct pair_class *pair2=pairs+b;
      if(pair1->i1->i==pair2->i1->i){
	i=pair1->i2->i; j=pair2->i2->i;
      }else if(pair1->i1->i==pair2->i2->i){
	i=pair1->i2->i; j=pair2->i1->i;
      }else if(pair1->i2->i==pair2->i1->i){
	i=pair1->i1->i; j=pair2->i2->i;
      }else if(pair1->i2->i==pair2->i2->i){
	i=pair1->i1->i; j=pair2->i1->i;
      }else{
	continue;
      }
      if(label_ij[i][j]<0)continue;
      pair=pairs+label_ij[i][j];
      //float P=pair1->score.P_C_nat*pair2->score.P_C_nat;
      float P=pair1->mut_inf*pair2->mut_inf;
      if((pair->ik<0)||(P > pair->P_ind)){
	pair->P_ind=P;
	pair->ik=a;
	pair->jk=b;
      }   
    }
  }
  pair=pairs;
  for(a=0; a<Nij_class; a++){
    if(pair->ik>=0){
      pair->P_ind=
	(pairs+pair->ik)->score.P_C_nat*(pairs+pair->jk)->score.P_C_nat;
    }
    pair++;
    }
}

void Print_CM(struct protein *prot, FILE *file_out, int *ali){
  fprintf(file_out, "# %d %d %s\n", prot->nres, prot->Nc, prot->name);
  for(int k=0; k<prot->Nc; k++){
    fprintf(file_out, "%d %d", prot->cont_list[k].res1,
	    prot->cont_list[k].res2);
    if(ali)fprintf(file_out, "   %d %d", ali[prot->cont_list[k].res1],
	    ali[prot->cont_list[k].res2]);
    fprintf(file_out, "\n");
  }
}

void Print_alignments(char *file_ali, int **ali_seq, char **MSA,
		      char **name, int *aligned, struct protein *template,
		      int N_template)
{
  char nameout[300]; Get_name(nameout, file_ali);
  strcat(nameout, ".ali");
  FILE *file_out=fopen(nameout, "w"); int k, l=0, line=60;
  printf("Writing alignments in file %s\n", nameout);
  for(int i=0; i<N_template; i++){
    int *ali=ali_seq[i], s=aligned[i];
    if(s<0)continue;
    struct protein *prot=template+i;
    fprintf(file_out, ">%s\n", prot->name); l=0;
    for(k=0; k<prot->nres; k++){
      fprintf(file_out, "%c", prot->amm[k]);
      l++; if(l==line){fprintf(file_out, "\n"); l=0;}
    }
    if(l)fprintf(file_out, "\n");
    fprintf(file_out, ">%s\n", name[s]); l=0;
    for(k=0; k<prot->nres; k++){
      if(ali[k]>=0){
	fprintf(file_out, "%c", MSA[s][ali[k]]); 
      }else{
	fprintf(file_out, "-");
      }
      l++; if(l==line){fprintf(file_out, "\n"); l=0;}
    }
    if(l)fprintf(file_out, "\n");
  }
  fclose(file_out);
  //exit(8);
}

int Seq_length(char *Seq, int L_ali){
  char *s=Seq; int l=L_ali;
  for(int j=0; j<L_ali; j++){if(*s=='-')l--; s++;}
  return(l);
}

void Copy_seq(char *Seq_k, char *Seq_i, int L_ali){
  char *sk=Seq_k, *si=Seq_i;
  for(int j=0; j<L_ali; j++){*sk=*si; si++; sk++;}
}

float Seq_id(char *Seq_i, char *Seq_j, int L_ali){
  int id=0; char *si=Seq_i, *sj=Seq_j;
  for(int j=0; j<L_ali; j++){
    if((*sj==*si)&&(*si!='-')){id++;} si++; sj++;
  }
  return((float)id/L_ali);
}

int Mut_min(char *Seq_i, char *Seq_j, int L_ali, int mutmin){
  // If mut > mutmin returns 1, otherwise returns 0 
  int m=0; char *si=Seq_i, *sj=Seq_j;
  for(int j=0; j<L_ali; j++){
    if((*sj!=*si)&&(*si!='-')&&(*sj!='-')){m++; if(m>=mutmin)return(1);}
    si++; sj++;
  }
  return(0);
}
