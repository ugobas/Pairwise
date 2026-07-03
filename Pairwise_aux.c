#include "pairwise.h"
#include "allocate.h"
#include "protein4.h"
#include "energy_BKV.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "Pairwise_aux.h"
#include "optimization.h"
#include "REM_new.h"
#include "externals.h"

// Parameters for selecting compact structures
#define C_OFFSET 4.0  // Asymptotic value of Nc/N
#define C_SLOPE 8.0   // Nc/N ~ C_OFFSET - C_SLOPE N^(-1/3)
//#define C_SLOPE 8.2   // Nc/N ~ C_OFFSET - C_SLOPE N^(-1/3)
#define C_THR 1.0      // Maximum deviation of Nc/N from typical value
                       // Very important parameter!!!
int L_MAX;
float ERR_THR=0.04; // Mean accepted error for Normalize_Q ERR_THR/40=0.001

// Energy parameters
int CP_Error=0;

// Number of contacts
extern int C_low, C_hig; 

// Concact statistics
float *Cont_f, *Cz1, *Cz2, *Cnc1;
float *nc_nc;
// <Nc>, <Nc^2>, <Nc^3>
//float *Nc1L, *Nc2L, *Nc3L;
float *nc2L;
float *Nc2L_center, *Nc3L_center;
double y1_ave, z_expo, Slope_x;

void Normalize_cont_freq(float *Cont_norm, float *Cont_f,
			 float *Cnc1_norm, float *Cnc1,
			 float *nc_nc_norm, float *nc_nc,
			 float *CNc1_norm, float *Cz1,
			 float *CNc2_norm, float *Cz2,
			 float *Cont_indir_norm, float *Cont_f_indir,
			 double Nc1L, double Nc2L, double Nc3L,
			 double NcL_indir, int L);

/*float Compute_cont_freq(int ij, int L);
float Compute_cont_freq_Nc(int ij, int L);
float Compute_cont_freq_Nc2(int ij, int L);
float Compute_cont_freq_indir(int ij, int L);*/

// Pairwise probabilities Q, predicted
float qmin=0.001, lmin; // Minimum value of q

void Compute_Q(float *Q, float Lambda, struct pair_class *pair, int cont);
float Normalize_Qprime_old(float *Q_pred, struct pair_class  *pair);
float Normalize_Qprime(float *Q_pred, float *Qprime,
		       struct pair_class  *pair, int cont);
float Normalize_Qprime_symm(float *Q_pred, float *Qprime,
			    struct pair_class  *pair, int cont);
void Update_q(float *qnew, float *Q, float *zi, float *zj);
float Update_z_naive(float *zi, float *zj, float *Q, struct pair_class *pair);
float Update_z_Newton(float *zi, float *zj, float *yi, float *yj,
		      float *Q, struct pair_class *pair);
float Error_norm(float *Q, struct pair_class *pair);
void Print_error_Q(struct pair_class *pair, float *Q,
		   float * zi, float * zj, float * ziopt, float * zjopt);
int printconv=0;


// Q observed
double Marginalize_N2_symm(float *P1, double *N2_obs);
void Pairwise_Q(float *Q_obs, float *P1i1, float *P1i2, double *N2_obs);

// Scores
float Contact_score(struct pair_class *pair, float *E_cont_not, int Na2);
void Compute_logarithm(float *log_Q, float *Q, float qmin, float lmin, int Na2);
float Compute_specific_mut_inf(struct pair_class *pair);
float dKL(float *log_Qpred, struct pair_class *pair);
float MSE(float *log_Q_pred, struct pair_class *pair);
float Likelihood(float *log_Qpred, struct pair_class *pair);
static float Corr_coeff(float *offset, float *slope, float *x, float *y, int n);
float Corr_coeff_w(float *offset, float *slope,
		   float *x, float *y, double *w, int n);

// Indirect contacts
float *Cont_freq_indir;
float NcL0_indir, NcL1_indir;    // <Nc_indirect>~L(NCL0+NCL1*L^(-1/3))
float *NcL_indir;

int Indirect_log_Q_PDB(float *log_Q, 
		       struct pair_class *pair_ij,
		       struct pair_class *pairs,
		       int Nij_class); //, int N_C
int Indirect_Q_new(float *Q, struct pair_class *pair_ij,
		   struct pair_class *pair_1, struct pair_class *pair_2,
		   float d_thr);
void Indirect_Q_ali(float *Q, int i_ij, struct pair_class *pairs);
void Indirect_Q_k(double *Q, int Naa,
		  struct pair_class *pair_ij,
		  struct pair_class *pair_ik,
		  struct pair_class *pair_jk);

/************************** General ***********************************/
float *Vectorize(float **V2, int Naa){
  float *V1=malloc((Naa*Naa)*sizeof(float)), *V=V1;
  for(int a=0; a<Naa; a++)for(int b=0; b<Naa; b++){*V=V2[a][b]; V++;}
  return(V1);
}

void Zero_mean(float *E_not, float *E,
	       float *P1i1, //pair->P1i1_obs
	       float *P1i2, //pair->P1i2_obs    struct pair_class *pair,
	       int Naa)
{
  float E1[Naa], E2[Naa], E12=0; int a,b;
  for(a=0; a<Naa; a++){E1[a]=0; E2[a]=0;}
  float *X=E;
  for(a=0; a<Naa; a++){
    for(b=0; b<Naa; b++){
      E1[a]+=P1i2[b]*(*X);
      E2[b]+=P1i1[a]*(*X);
      X++;
    }
    E12+=P1i1[a]*E1[a];
  }
  float *Y=E_not; X=E;
  for(a=0; a<Naa; a++){
    for(b=0; b<Naa; b++){
      *Y=*X-E1[a]-E2[b]+E12; X++; Y++;
    }
  }
} 


/************************ Contact frequencies ***************************/

// Contact statistics

void Initialize_REM_2(int L_max, int ij_min,
		      struct protein *prot_ptr, int N_prot)
{
  IJ_MIN=ij_min;
  L_MAX=L_max;

  // Allocate
  Cont_f =malloc(L_max*sizeof(float));
  Cnc1  =malloc(L_max*sizeof(float));
  Cz1  =malloc(L_max*sizeof(float));
  Cz2  =malloc(L_max*sizeof(float));
  nc_nc =malloc(L_max*sizeof(float));

  printf("### Contact statistics\n%d contact matrices read\n", N_prot);

  /* Scaling of contacts */
  float Slope_tmp=C_SLOPE/C_OFFSET; 
  float Ncp[N_prot], Xp[N_prot], logLp[N_prot];
  int np=0, select[N_prot], i;
  for(i=0; i< N_prot; i++){
    struct protein *prot=prot_ptr+i;
    float L=prot->nres;
    float x=pow(L, -0.333333), c=(prot->Nc/L);
    if(fabs(c-C_OFFSET*(1-Slope_tmp*x))>C_THR){
      select[i]=0;
    }else{
      select[i]= 1;
      logLp[np]=log(L); Xp[np]=x; Ncp[np]=c;
      np++;
    }
  }
  // Fits of Nc(L), Nc^2(L) etc.
  printf("Fitting number of contacts versus surface to volume scaling.\n");
  printf("%d proteins kept and %d omitted due to large deviations\n",
	 np, N_prot-np);
  // First moment
  float Nc1L0, Nc1L1;    // <Nc>~L(NC1L0+NC1L1*L^(-1/3))
  float r=Corr_coeff(&Nc1L0, &Nc1L1, Xp, Ncp, np);
  printf("<Nc/L> ~ %.2f + %.1fL^(-1/3)  r= %.3f\n", Nc1L0, Nc1L1, r);
  Slope_x=Nc1L1/Nc1L0;

  // Normalize Nc/L by fitted value
  float Yp[np]; y1_ave=0;
  printf("y = NC/L(1 + %.2fL^(-1/3))\n", Slope_x);
  for(i=0; i< np; i++){
    Yp[i]=Ncp[i]/(1+Slope_x*Xp[i]);
    y1_ave+=Yp[i];
  } 
  y1_ave/=np;
  float coeff, slope;
  r=Corr_coeff(&coeff, &slope, Xp, Yp, np);
  printf("y ~ %.3f + %.1f*L^(-1/3)  r= %.3f\n", coeff, slope, r);
  printf("<y>= %.3f\n", y1_ave);

  // Second and third moment with power law
  float logY2p[np], logY3p[np];
  for(i=0; i< np; i++){
    float ly=log(fabs(Yp[i]-y1_ave));
    logY2p[i]=2*ly; logY3p[i]=3*ly;
  }   
  float Y2L0, Y2L1, r2;  // <(y-<y>)^2>~exp(Y2L0+Y2L1*log(L))
  r2=Corr_coeff(&Y2L0, &Y2L1, logLp, logY2p, np);
  printf("<(y-<y>)^2> ~ exp(%.2f + %.1flog(L))  r= %.3f\n", Y2L0,Y2L1,r2);
  float Y3L0, Y3L1, r3;  // <(y-<y>)^3>~exp(Y3L0+Y3L1*log(L))
  r3=Corr_coeff(&Y3L0, &Y3L1, logLp, logY3p, np);
  printf("<(y-<y>)^3> ~ exp(%.2f + %.1flog(L))  r= %.3f\n", Y3L0,Y3L1,r3);
  r2*=r2; r3*=r3;
  z_expo=(r2*Y2L1/2+r3*Y3L1/3)/(r2+r3);
  printf("Average scaling exponent: %.3f\n", z_expo);

  // Compute Z score
  float Zp[np]; double z2_ave=0, z3_ave=0; 
  for(i=0; i< np; i++){
    float z=(Yp[i]-y1_ave)*exp(-z_expo*logLp[i]);
    Zp[i]=z; z2_ave+=z*z; z3_ave+=z*z*z;
  }
  z2_ave/=np; z3_ave/=np;
  r=Corr_coeff(&coeff, &slope, Xp, Zp, np);
  printf("z=(y-<y>)*L^%.2f ~ %.2f + %.2fL^(-1/3)  r= %.3f\n",
	 -z_expo,coeff,slope,r);
  printf("<z^2>= %.3g\n", z2_ave);

  /**************** Statistics of contacts ****************/
  long Cont_num[L_max], Cont_norm[L_max], nc_nc_num[L_max],
    ncCont_num[L_max], z_Cont_num[L_max], z2_Cont_num[L_max];
  for(i=0; i<L_max; i++){
    Cont_num[i]=0; Cont_norm[i]=0;
    ncCont_num[i]=0; nc_nc_num[i]=0;
    z_Cont_num[i]=0; z2_Cont_num[i]=0;
  }

  float nc2[N_prot]; int kp=0;
  for(int ip=0; ip< N_prot; ip++){
    if(select[ip]==0)continue;
    struct protein *prot=prot_ptr+ip;
    int L=prot->nres;
    int nc[L], ncont_l[L];
    for(i=0; i<L; i++){nc[i]=0; ncont_l[i]=0;}
    for(int ic=0; ic<prot->Nc; ic++){
      struct contact *cont=prot->cont_list+ic;
      ncont_l[cont->res2-cont->res1]++;
      nc[cont->res1]++; nc[cont->res2]++;
    }
    // <c_ij ni>
    for(int ic=0; ic<prot->Nc; ic++){
      struct contact *cont=prot->cont_list+ic;
      ncCont_num[cont->res2-cont->res1]+=(nc[cont->res1]+nc[cont->res2]);
    }
    // Sum protein to global counters
    float z=Zp[kp];
    float scale2=(1+Slope_x*Xp[kp]); scale2*=scale2;
    //float z=Yp[kp];
    for(i=0; i<L; i++){
      if(i>=L_max)break;
      Cont_norm[i]+=L-i;
      Cont_num[i]+=ncont_l[i];
      z_Cont_num[i]+=ncont_l[i]*z;  // Average contacts per residue
      z2_Cont_num[i]+=ncont_l[i]*z*z;
      for(int j=0; j<(L-i); j++){
	if(j >= L_max)break;
	nc_nc_num[i]+=nc[j]*nc[j+i]/scale2;
      }
    } // End pairs
    // sum_i (ni/L)^2
    float n2=0; for(i=0; i<L; i++)n2+=nc[i]*nc[i]/scale2;
    nc2[kp]=n2/L; kp++;
  } // end proteins
  if(kp!=np){
    printf("ERROR, wrong number of proteins %d instead of %d\n",
	   kp, np); exit(8);
  }

  float nc2L0, nc2L1;  // sum_i <(ni/L)^2>~(nc2L0+nc2L1*L^(-1/3))
  r=Corr_coeff(&nc2L0, &nc2L1, Xp, nc2, np);
  printf("sum_i ni^2/L ~ %.2f + %.1fL^(-1/3)  r= %.3f\n", nc2L0, nc2L1, r);
  nc2L0=0; for(i=0; i<np; i++)nc2L0+=nc2[i]; nc2L0/=np;
  printf("Average sum_i ni^2/L: %.2f\n", nc2L0);

  //printf("Cont_freq <NcCij>\n");
  for(i=0; i<L_max; i++){
    float norm=Cont_norm[i];
    if(Cont_num[i]>norm)norm=Cont_num[i];
    Cont_f[i]=(float)Cont_num[i];
    Cnc1[i]=(float)ncCont_num[i];
    Cz1[i]=(float)z_Cont_num[i];
    Cz2[i]=(float)z2_Cont_num[i];
    nc_nc[i]=nc_nc_num[i];

    if(norm){
      Cont_f[i]/=norm;
      Cnc1[i]/=(2*norm);
      Cz1[i]/=norm;
      Cz2[i]/=norm;
      nc_nc[i]/=norm;
    }else{
      printf("WARNING, i=%d norm=0\n", i);
    }

    //if((i>=4)&&(i<30))printf("%.4f %.3f\n", Cont_f[i], CNc1[i]);
    if(Cont_f[i]<0){
      printf("ERROR, negative Cont_freq:\n");
      printf("%d cf= %.4f %ld %ld Ncf= %.4f Nc2f= %.4f\n",
	     i, Cont_f[i], Cont_num[i], Cont_norm[i], Cz1[i], Cz2[i]);
      exit(8);
    }
  }

  nc2L=malloc(L_max*sizeof(float));
  Nc1L=malloc(L_max*sizeof(float));
  Nc2L=malloc(L_max*sizeof(float));
  Nc3L=malloc(L_max*sizeof(float));
  Nc2L_center=malloc(L_max*sizeof(float));
  Nc3L_center=malloc(L_max*sizeof(float));
  //printf("#L N<Nc> <Nc^2> <Nc^2>-<Nc>^2\n");
  for(i=4; i<L_max; i++){
    double x=pow(i,-0.33333333);
    float core=(1+Slope_x*x);
    //nc2L[i]=(nc2L0+nc2L1*x)*i;
    nc2L[i]=(nc2L0*core*core)*i;
    double scale=pow(i, z_expo);
    float Lcore=core*i;
    Nc1L[i]=y1_ave*Lcore;
    float y2=y1_ave*y1_ave, z2=z2_ave*scale*scale;
    Nc2L[i]=(y2+z2)*Lcore*Lcore;
    float Lscale=Lcore*scale;
    Nc2L_center[i]=z2_ave*Lscale*Lscale;
    Nc3L_center[i]=z3_ave*Lscale*Lscale*Lscale;
    Nc3L[i]=y1_ave*(y2+3*z2)*Lcore*Lcore*Lcore+Nc3L_center[i];
    //if((i>90)&&(i<120))
    //printf("%d %.3g %.3g %.3g\n", i, Nc1L[i], Nc2L[i], Nc2L_center[i]);
  }
  //exit(8);

}

void Indirect_cont_stat(int L_max, int ij_min,
			struct protein *prot_ptr, int N_prot)
{
  // Allocate
  int i;
  Cont_freq_indir =malloc(L_max*sizeof(float));
  long *Cont_num=malloc(L_max*sizeof(long));
  long *Cont_norm=malloc(L_max*sizeof(long));
  for(i=0; i<L_max; i++){
    Cont_num[i]=0; Cont_norm[i]=0;
  }

  /* Statistics of contacts */
  float *Nc1=malloc(N_prot*sizeof(float));
  float *Nc2=malloc(N_prot*sizeof(float));
  float *Xp=malloc(N_prot*sizeof(float));
  for(int ip=0; ip< N_prot; ip++){
    struct protein *prot=prot_ptr+ip;
    float c=prot->Nc_indirect/(float)prot->nres;
    Nc1[ip]=c; Nc2[ip]=c*c;
    Xp[ip]=pow(prot->nres, -0.333333);
    if(prot->indirect_contacts==NULL)continue;
    struct contact *cont=prot->indirect_contacts;
    for(i=0; i<prot->Nc; i++){
      int l=cont->res2-cont->res1; cont++;
      if(l<L_max)Cont_num[l]++; 
    }
    for(i=0; i<L_max; i++){
      int nl=prot->nres-i;
      if(nl<=0)break;
      Cont_norm[i]+=nl;
    }
  }
  printf("Cont_freq_indirect:\n");
  for(i=0; i<L_max; i++){
    Cont_freq_indir[i]=(float)Cont_num[i]/(float)Cont_norm[i];
    if(i && i<40)printf("%d %.4f\n", i, Cont_freq_indir[i]);
    if((Cont_freq_indir[i]<0)){
      printf("ERROR, negative Cont_freq:\n");
      printf("%d cf= %.4f %ld %ld \n",
	     i, Cont_freq_indir[i], Cont_num[i], Cont_norm[i]);
      exit(8);
    }
  }
  float r=Corr_coeff(&NcL0_indir, &NcL1_indir, Xp, Nc1, N_prot);
  printf("NC_indirect/L ~ %.2f + %.2f^(-1/3)  r= %.3f\n",
	 NcL0_indir, NcL1_indir, r);

  NcL_indir=malloc(L_max*sizeof(float));
  for(i=10; i<L_max; i++){
    float x=pow(i,-0.33333);
    NcL_indir[i]=(NcL0_indir+NcL1_indir*x)*i;
  }
}


/*
float Compute_cont_freq(int ij, int L){
  double norm=0;
  for(int l=IJ_MIN; l<L; l++)norm+=Cont_f[l]*(L-l);
  return(Cont_f[ij]*Nc1L[L]/norm);
}

float Compute_cont_freq_Nc(int ij, int L){
  double norm=0;
  for(int l=IJ_MIN; l<L; l++)norm+=CNc1[l]*(L-l);
  return(CNc1[ij]*Nc2L[L]/norm);
}

float Compute_cont_freq_Nc2(int ij, int L){
  double norm=0;
  for(int l=IJ_MIN; l<L; l++)norm+=CNc2[l]*(L-l);
  return(CNc2[ij]*Nc3L[L]/norm);
}

float Compute_cont_freq_indir(int ij, int L){
  double norm=0;
  for(int l=IJ_MIN; l<L; l++)norm+=Cont_freq_indir[l]*(L-l);
  return(Cont_freq_indir[ij]*NcL_indir[L]/norm);
}
*/

void Normalize_cont_freq(float *Cont_norm, float *Cont_f,
			 float *Cnc1_norm, float *Cnc1,
			 float *nc_nc_norm, float *nc_nc,
			 float *CNc1_norm, float *Cz1,
			 float *CNc2_norm, float *Cz2,
			 float *Cont_indir_norm, float *Cont_f_indir,
			 double Nc1, double Nc2, double Nc3,
			 double Nc_indir, int L)
{
  double sum_c=0, sum_c_nc=0, sum_c_z=0,
    sum_c_z2=0, sum_nc_nc=0, sum_nl=0, sum_c_indir=0;
  int l;
  for(l=IJ_MIN; l<L; l++){
    int nl=L-l;
    sum_nl+=nl;
    sum_c+=Cont_f[l]*nl;
    sum_c_nc+=Cnc1[l]*nl;
    sum_c_z+=Cz1[l]*nl;
    sum_c_z2+=Cz2[l]*nl;
    sum_nc_nc+=nc_nc[l]*nl;
    if(Cont_f_indir)sum_c_indir+=Cont_f_indir[l]*nl;
  }
  if(sum_c==0){
    printf("ERROR in normalization sum_c, nan\n"); exit(8);
  }
  if(sum_c_z==0){
    printf("ERROR in normalization sum_c_z, nan\n"); exit(8);
  }
  if(sum_c_z2==0){
    printf("ERROR in normalization sum_c_z2, nan\n"); exit(8);
  }
  if(sum_c_nc==0){
    printf("ERROR in normalization sum_c_nc, nan\n"); exit(8);
  }
  if(sum_nc_nc==0){
    printf("ERROR in normalization sum_nc_nc, nan\n"); exit(8);
  }

  float Nc12=Nc1*Nc1;
  float norm_c_z= (Nc2-Nc12)/sum_c_z;
  float norm_c_z2=(Nc3-Nc12*Nc1)/sum_c_z2;
  float norm_c=Nc1/sum_c;
  float norm_c_nc=nc2L[L]/sum_c_nc;
  float norm_nc_nc=2*Nc2/sum_nc_nc; // sum_i<j n_i*nj=4Nc/2!!
  float norm_c_indir=0;
  if(Cont_f_indir)norm_c_indir=Nc_indir/sum_c_indir;
  float nc=2*Nc1/L, nc2=2*Nc1*Nc1/sum_nl;

  // float cc=Nc1L*Nc1L/norm_nl;
  for(l=IJ_MIN; l<L; l++){
    Cont_norm[l]=Cont_f[l]*norm_c;
    if(Cont_norm[l]>1)Cont_norm[l]=1;
    nc_nc_norm[l]=nc_nc[l]*norm_nc_nc-nc2;
    CNc1_norm[l]=Cz1[l]*norm_c_z;
    CNc2_norm[l]=Cz2[l]*norm_c_z2-2*Nc1*CNc1_norm[l];
    Cnc1_norm[l]=Cnc1[l]*norm_c_nc-Cont_norm[l]*nc;
    if(Cont_f_indir)Cont_indir_norm[l]=Cont_f_indir[l]*norm_c_indir;
  }
}

void  Initialize_prot_class(struct prot_class *prot){
  prot->Np=0;
  prot->L=0;
  prot->Nc=0;
  prot->npairs=0;
  prot->U_ave=0;
  //for(int a=0; a<Naa; a++)prot->U1[a]=0;
  prot->U_norm=0;
  for(int a=0; a<Naa; a++)prot->P1_glob[a]=0;
  prot->Q_global=malloc(Na2*sizeof(float));
  prot->Q_glob_ini=malloc(Na2*sizeof(float));
  prot->Q_glob_opt=malloc(Na2*sizeof(float));
  prot->log_Q_glob=malloc(Na2*sizeof(float));
  for(int c=0; c<3; c++)prot->Npair[c]=0;
}

void Initialize_pair(struct pair_class *pair, struct prot_class *prot)
{
  pair->prot_class=prot;
  // Set to zero
  pair->ij=0;
  pair->Cont_freq=0;
  pair->nc_nc=0;
  pair->Cont_freq_nc=0;
  pair->Cont_freq_Nc=0;
  pair->Cont_freq_Nc2=0;
  pair->Cont_freq_indir=0;
  pair->E_cont_obs=0;
  pair->min_dist=-1;
  pair->n=0;

  pair->N2_obs=malloc(Na2*sizeof(double));
  for(int a=0; a<Na2; a++)pair->N2_obs[a]=0;
  pair->Q_obs=malloc(Na2*sizeof(float));
  pair->log_Q_obs=malloc(Na2*sizeof(float));
  for(int i=0; i<3; i++){
    pair->log_Q_pred[i]=malloc(Na2*sizeof(float));
    pair->Q_pred[i]=malloc(Na2*sizeof(float));
  }

  pair->ene=malloc(2*sizeof(float *));
  for(int c=0; c<=1; c++){
    pair->ene[c]=malloc(Na2*sizeof(float));
  }

}

void Empty_pair(struct pair_class *pair)
{
  free(pair->N2_obs);
  free(pair->Q_obs);
  free(pair->log_Q_obs);
  for(int c=0; c<3; c++){
    free(pair->log_Q_pred[c]);
    free(pair->Q_pred[c]);
  }
}


void Initialize_site(struct site *site, struct prot_class *prot, int i)
{
  site->prot_class=prot;
  site->i=i;
  site->n=0;
  site->nc=0;
  site->i_c=i;
  for(int a=0; a<21; a++){
    site->P1_obs[a]=0;
    site->U1[a]=0;
    site->cU1[a]=0;
  }
}

int Set_index_nc(int nc, int N_C)
{
  int i_c;
  if(nc<=C_low){return(0);}
  else if(nc>C_hig){i_c=2;}
  else{return(1);}
  if(i_c>=N_C){i_c=N_C-1;}
  return(i_c);
}

int Set_index_nc2(int i1, int i2, int N_C)
{
  if(i1==0){return(i2);}
  else if(i1==1){return(2+i2);}
  else{return(5);}
  // int ncp=(N_C*i1)-i1+(i2-i1)=(N_C-2)*i1+i2 
}

int Sum_pairs(struct prot_class *prot_c, int Npair_prot,
	      struct protein *protein, int L_seq,
	      int **C_mat, int **C2, int N_ij, int *ij_bin,
	      int ij_min, int N_C, int **label_ij)
{
  //int N_C2=N_C*(N_C+1)/2; // Number of pairs of sites (6 for N_C=3)
  int Npairs=0, L=protein->nres;
  prot_c->L+=L; 
  prot_c->Nc+=Nc1L[L];
  float wp=protein->Nc;
  prot_c->U_ave+=wp*protein->U_ave;
  prot_c->U_norm+=wp;
  prot_c->Np++;

  // Normalization of contacts
  float Cont_norm[L], Cnc1_norm[L], CNc1_norm[L], CNc2_norm[L],
    Cont_indir_norm[L], nc_nc_norm[L];
  Normalize_cont_freq(Cont_norm, Cont_f,
		      Cnc1_norm, Cnc1,
		      nc_nc_norm, nc_nc,
		      CNc1_norm, Cz1,
		      CNc2_norm, Cz2,
		      Cont_indir_norm, Cont_freq_indir,
		      Nc1L[L], Nc2L[L], Nc3L[L], NcL_indir[L], L);

  int k=0, *nc=protein->n_cont;
  int *i_aa=protein->i_aa;
  struct pair_class *pair_c=prot_c->pairs, *pair=pair_c;
  for(int i=0; i<L_seq; i++){

    int i_c=-1;
    if(label_ij==NULL){i_c=Set_index_nc(nc[i], N_C);} // PDB struct
    struct site *site=prot_c->site+i_c;
    int a=i_aa[i];
    if(a<0){if(Naa<=20){continue;} a=20;}
    if(a<=20)site->P1_obs[a]++;
    int na=a*Naa;

    site->n++;
    site->nc+=nc[i];
 
    for(int j=i+ij_min; j<L_seq; j++){
      int b=i_aa[j]; if(b<0){if(Naa<=20){continue;} b=20;}
      int ij=j-i; if(ij>=L)ij=L-1;
      int C;
      if(C_mat){C=C_mat[i][j]; if(C2[i][j] && C==0){C=2;}}
      else{C=-1;}
      if(label_ij==NULL){
	int i1, i2, j_c=Set_index_nc(nc[j], N_C);
	if(i_c<=j_c){i1=i_c; i2=j_c;}else{i1=j_c; i2=i_c;}
	// PDB structures: C_nat=0, 1, 2
	/*k=Find_ij_class(ij, C_mat[i][j], C2[i][j], N_ij, ij_bin,
	  i1, i2, N_C, N_C2);
	  pair=pair_c+k; */
	pair=Find_pair(pair_c, Npair_prot, i1, i2, C, ij, N_ij, ij_bin);
	pair->i1=prot_c->site+i1;
	pair->i2=prot_c->site+i2;

	if(i_c==j_c){
	  // Symmetrize observed distribution
	  pair->N2_obs[na+b]+=0.5;
	  pair->N2_obs[Naa*b+a]+=0.5;
	}else if(i_c<j_c){
	  pair->N2_obs[na+b]++;
	}else{
	  pair->N2_obs[Naa*b+a]++;
	}

      }else { //if(label_ij){
	// MSA

	// Remind to eliminate columns not present in target
	pair=pair_c+k; 
	pair->N2_obs[na+b]++;
	pair->i1=prot_c->site+i;
	pair->i2=prot_c->site+j;
	pair->ij_index=j-i;
	pair->C_nat=C;
	k++;
      }

      Npairs++;
      pair->n++;
      pair->ij+=ij;
      pair->E_cont_obs+=E_cont_T[a][b];

      pair->Cont_freq+=Cont_norm[ij];
      pair->Cont_freq_nc+=Cnc1_norm[ij];
      pair->Cont_freq_Nc+=CNc1_norm[ij];
      pair->Cont_freq_Nc2+=CNc2_norm[ij];
      pair->nc_nc+=nc_nc_norm[ij];
      pair->Cont_freq_indir+=Cont_indir_norm[ij];
    }
  }
  prot_c->npairs+=Npairs;
  return(Npairs);
}

void Normalize_sites(struct prot_class *prot, int Nsites_prot, int Npair_prot)
{
  struct site *sites=prot->site; int a, b;
  for(int k=0; k<Nsites_prot; k++){
    struct site *site=sites+k;
    site->nc/=site->n;
    double sum=0;
    for(a=0; a<Naa; a++){sum+=site->P1_obs[a];}
    for(a=0; a<Naa; a++){site->P1_obs[a]/=sum;}
    for(a=0; a<Naa; a++){
      site->U1[a]=0;
      for(b=0; b<Naa; b++){
	site->U1[a]+=E_cont_T[a][b]*site->P1_obs[b];
      }
    }
  }
  for(int k=0; k<Nsites_prot; k++){
    struct site *site=sites+k;
    for(a=0; a<Naa; a++){site->cU1[a]=0;}
    for(int l=0; l<Npair_prot; l++){
      struct pair_class *p= prot->pairs+l;
      if(p->i1 == site){
	for(a=0; a<Naa; a++){
	  site->cU1[a]+= p->Cont_freq* (p->i2)->U1[a];
	}
      }
      if(p->i2 == site){
	for(a=0; a<Naa; a++){
	  site->cU1[a]+= p->Cont_freq* (p->i1)->U1[a];
	}
      }
    }
  }
}

void Protein_stat(struct prot_class *prot, int Npair_prot,
		  int ij_min, int L_max, int I_WEIGHT)
{ 
  int TEST=1;
  if(prot->Np==0)return;
  if(prot->Np>1){
    prot->L/=prot->Np;
    prot->Nc/=prot->Np;
  }
  prot->U_ave/=prot->U_norm;
  //for(a=0; a<Naa; a++)prot->U1[a]/=prot->U_norm;
  //if(Naa==21)prot->U1[20]=E_cont_T[20][20];

  double norm_c=0, norm_nc_nc=0, norm_c_nc=0, //norm=0, 
    norm_c_Nc=0, norm_c_Nc2=0, norm_nl=0;
  struct pair_class *pair; int j;
  for(j=0; j<Npair_prot; j++){
    pair=prot->pairs+j;
    if(pair->n<=NTHR)continue;
    prot->Npair[pair->C_nat]++;
    if(pair->C_nat==1){pair->min_dist=4.4;}
    else if(pair->C_nat==2){pair->min_dist=6.0;}
    else{pair->min_dist=15;}
    pair->mut_inf=Compute_specific_mut_inf(pair);
    if(I_WEIGHT==0){pair->w=1;}
    else if(I_WEIGHT==1){pair->w=pair->n;}
    else if(I_WEIGHT==2){pair->w=log(pair->n);}
    else if(I_WEIGHT==3){pair->w=pair->mut_inf;}
    //norm+=pair->w;

    pair->E_cont_obs/=pair->n;
    pair->ij/=pair->n;
    pair->Cont_freq/= pair->n;
    pair->Cont_freq_nc/= pair->n;
    pair->Cont_freq_Nc/= pair->n;
    pair->Cont_freq_Nc2/= pair->n;
    pair->Cont_freq_indir/=pair->n;
    pair->nc_nc/= pair->n;

    if(pair->ij < ij_min){
      printf("ERROR, too small <|i-j|>= %.1f < %d n=%.0f\n",
	     pair->ij, ij_min, pair->n);
      exit(8);
    }else if(pair->ij >= L_max){
      printf("Protein stat |i-j|= %.1f -> %d\n", pair->ij, L_max-1);
      pair->ij=L_max-1;
    }
    if(prot->Np){
      float nl= pair->n/prot->Np; // Average number of pairs per protein
      pair->num_cont=nl;
      norm_nl+=nl;
      norm_c+=pair->Cont_freq*nl;
      norm_nc_nc+=pair->nc_nc*nl;
      norm_c_nc+=pair->Cont_freq_nc*nl;
      norm_c_Nc+=pair->Cont_freq_Nc*nl;
      norm_c_Nc2+=pair->Cont_freq_Nc2;
    }
  }

  if(TEST){
    printf("Prot: L=%.1f Nc=%.f n=%d\n", prot->L, prot->Nc, prot->Np);
    printf("Prior to normalization:\n");
    for(j=0; j<Npair_prot; j++){
      pair=prot->pairs+j;
      printf("pair %d i=%d j=%d C=%d ij=%.3g <C>=%6.3g"
	     " nc_nc=%6.3g C_NC= %.3f C_Nc2= %.4g"
	     "\tn=%8.0f\n",
	     j, pair->i1->i_c, pair->i2->i_c, pair->C_nat, pair->ij,
	     pair->Cont_freq, pair->nc_nc,
	     pair->Cont_freq_Nc, pair->Cont_freq_Nc2,
	     pair->n);
    }
  }

  pair=prot->pairs;
  int L=prot->L;
  norm_c=Nc1L[L]/norm_c;
  norm_nc_nc=2*Nc2L_center[L]/norm_nc_nc;
  norm_c_nc=0.5*(nc2L[L]-4*Nc1L[L]*Nc1L[L]/L)/norm_c_nc;
  norm_c_Nc=Nc2L_center[L]/norm_c_Nc;
  norm_c_Nc2=Nc3L_center[L]/norm_c_Nc2;
  for(j=0; j<Npair_prot; j++){
    pair=prot->pairs+j;
    if(pair->n<=NTHR)continue;
    pair->Cont_freq *= norm_c;
    pair->nc_nc *= norm_nc_nc;
    pair->Cont_freq_nc *=norm_c_nc;
    pair->Cont_freq_Nc *=norm_c_Nc;
    pair->Cont_freq_Nc2 *=norm_c_Nc2;
    //Compute_Phi(pair, prot);
    Compute_ene(pair, prot);
  }

  if(TEST){
    printf("After normalization:\n");
    for(j=0; j<Npair_prot; j++){
      pair=prot->pairs+j;
      printf("pair %d i=%d j=%d C=%d ij=%.3g <C>=%6.3g"
	     " nc_nc=%6.3g C_NC= %.3f C_Nc2= %.4g"
	     "\tn=%8.0f\n",
	     j, pair->i1->i_c, pair->i2->i_c, pair->C_nat, pair->ij,
	     pair->Cont_freq, pair->nc_nc,
	     pair->Cont_freq_Nc, pair->Cont_freq_Nc2,
	     pair->n);
    }
  }
}

/****************** Predicted Q *****************************/

void Compute_Qprime_old(float *Qprime, struct pair_class  *pair, int cont,
			float Lambda)
{
  // Wrong: Q' must not be symmetric. It must be "normalized" 
  /*printf("Computing beta, Phi1= %.3g Phi2= %.3g Lambda= %.3f\n",
    Phi1,Phi2,Lambda);*/
  float LPhi0=0, LPhi11=0, LPhi01=0, LPhi2=0;

  //if(REM){LPhi0=Lambda*(cont+pair->Phi0);}
  //LPhi11=Lambda*pair->Phi11;
  //LPhi01=Lambda*pair->Phi01, LPhi2=Lambda*pair->Phi2;
  float *U1=pair->i1->U1;

  for(int a=0; a<Naa; a++){
    int ab=a*Naa; float *q=Qprime+ab; 
    for(int b=0; b<=a; b++){
      double beta=E_cont_T[a][b]*LPhi0;
      if(REM>=2)
	beta+=E_cont_T[a][b]*(LPhi2*E_cont_T[a][b]+ LPhi01*(U1[a]+U1[b]))
	  +LPhi11*U1[a]*U1[b];
 
      *q=exp(-beta);
     if(isnan(*q) || isinf(*q)){
	printf("ERROR in compute_qprime, LPhi0= %.2g LPhi2= %.2g"
	       " LPhi01= %.2g LPhi11=%.2g beta=%.2g Lambda= %.2g\n",
	       LPhi0, LPhi2, LPhi01, LPhi11, beta, Lambda);
	exit(8);
      }
      if(b<a){int ba=b*Naa+a; Qprime[ba]=*q;} q++;
    }
  }
}

void Compute_ene(struct pair_class *pair, struct prot_class *prot)
{
  float Phi1=0, Phi2=0, Phi11=0;
  if(REM){Phi1=-pair->Cont_freq;}
  if(REM>=2){
    if(prot->U_ave>0){Phi1+=(pair->Cont_freq_Nc)*prot->U_ave;} //
    //pair->Phi0+=(pair->Cont_freq_Nc-pair->Cont_freq_nc)*prot->U_ave; //
    Phi11=0.25*pair->nc_nc/TEMP; 
    //float Cij_2=pair->Cont_freq*(1-pair->Cont_freq);
    //pair->Phi01=0.5*(pair->Cont_freq_nc-Cij_2)/TEMP;
    //Phi2=0.5*Cij_2/TEMP;
    Phi2=0.5*pair->Cont_freq*(1-pair->Cont_freq)/TEMP;
  }

  int error=0;
  if(isnan(Phi1)) {printf("ERROR, Phi1 is nan\n"); error++;}
  if(isnan(Phi11)){printf("ERROR, Phi11 is nan\n"); error++;}
  if(isnan(Phi2)){printf("ERROR, Phi2 is nan\n"); error++;}

  if(error || abs(Phi1) > 10 || abs(Phi11) >1000 || abs(Phi2)>1000){
    int L=prot->L;
    printf("WARNING, nan or very large Phi1= %.2f\n", Phi1);
    printf("L= %d <U>= %.2g ij= %.2g ", L, prot->U_ave, pair->ij);
    printf("cij= %.2f\n", pair->Cont_freq);
    printf("<Nc>= %.0f <Cij*Nc>-<Cij><Nc>= %.2f\n",
	   Nc1L[L], pair->Cont_freq_Nc);
    printf("C_ij_nc= %.3f\n", pair->Cont_freq_nc);
    printf("nc_nc= %.3f\n", pair->nc_nc);
    if(REM>=3)
      printf("<Nc^2>= %.0f Cij_Nc2= %.2g\n",
	       Nc2L[L], pair->Cont_freq_Nc2);
    exit(8);
  }

  int repel=0;
  if(REM>=2){
    float E_1=0, E_2=0;
    for(int a=0; a<Naa; a++){
      E_1+=pair->P1i1_obs[a]*pair->i1->cU1[a];
      E_2+=pair->P1i2_obs[a]*pair->i2->cU1[a];
    }
    if(E_1 >0 || E_2 >0){repel=1;}
  }
  int ab=0;
  for(int a=0; a<Naa; a++){
    for(int b=0; b<Naa; b++){ // ene[0]: no contact
      pair->ene[0][ab]=(Phi1+Phi2*E_cont_T[a][b])*E_cont_T[a][b];
      if(REM>=2 && repel==0){
	pair->ene[0][ab]+=Phi11*
	  (pair->i1->U1[a]*pair->i2->U1[b]+  // U1 or cU1?
	   pair->i1->U1[b]*pair->i2->U1[a]);
      }
      pair->ene[1][ab]=pair->ene[0][ab]+E_cont_T[a][b]; // ene[1]: contact
      ab++;
    }
  }  
}

void Compute_Q(float *Q, float Lambda, struct pair_class *pair, int cont)
{
  // log Q ~ log(Q_glob)+Lambda*ene_iajb
  float Ene_not[Na2];
  Zero_mean(Ene_not, pair->ene[cont], pair->P1i1_obs, pair->P1i2_obs, Naa);

  float *Q_glob=pair->prot_class->Q_global;

  float Qprime[Na2], *Qp=Qprime, *e_not=Ene_not;
  for(int a=0; a<Naa; a++){
    for(int b=0; b<Naa; b++){
      *Qp=exp(-Lambda*(*e_not));
      if(Q_GLOB){(*Qp) *= (*Q_glob); Q_glob++;}
      if(isnan(*Qp) || isinf(*Qp)){
	printf("ERROR in compute_Q, a=%d b=%d ene= %.2g q=%.2g Lambda= %.2g\n",
	       a, b, *e_not, *Qp, Lambda);
	exit(8);
      }
      Qp++; e_not++;
    }
  }
  Normalize_Qprime(Q, Qprime, pair, cont);
 
  return;
}

float Normalize_Qprime(float *Q_pred, float *Qprime,
		       struct pair_class  *pair, int cont)
{
  if(pair->i1==pair->i2){
    return(Normalize_Qprime_symm(Q_pred, Qprime, pair, cont));
  }
  
  int IT_MAX=500; float EPS=0.0001, EPS_MAX=0.01, e_min=1000;
  int a, b, iter;
  double z1, zij=0, zi[Naa], zj[Naa], zinew[Naa], zjnew[Naa];
  float ziopt[Naa], zjopt[Naa];
  // Inizialization
  for(a=0; a<Naa; a++){zi[a]=1; zj[a]=1;}
  float *q =Qprime, *P1i1=pair->P1i1_obs, *P1i2=pair->P1i2_obs;
  // Normalize Qprime
  for(a=0; a<Naa; a++){
    z1=0; float *p=P1i2;
    for(b=0; b<Naa; b++){z1+=(*p)*(*q); q++; p++;}
    zij+=P1i1[a]*z1;
  }
  zij=fabs(zij);
  if(zij<=0 || isnan(zij)){
    printf("ERROR in Normalize_Qp, zij= %.2g\n", zij);
    return(-1);
  }
  zij=1./zij; q=Qprime;
  for(a=0; a<Naa; a++){
    for(b=0; b<Naa; b++){(*q)*=zij; q++;}
  }
  
  float dsumi, dsumj;
  for(iter=0; iter<IT_MAX; iter++){
    // Compute zi, zj 
    dsumi=0; dsumj=0;
    float zPi[Naa], zPj[Naa];
    for(a=0; a<Naa; a++){
      zPi[a]=zi[a]*P1i1[a];
      zPj[a]=zj[a]*P1i2[a];
      zjnew[a]=0;
    }
    
    q=Qprime; int error=0;
    for(a=0; a<Naa; a++){
      z1=0; 
      for(b=0; b<Naa; b++){
	z1+=zPj[b]*(*q); zjnew[b]+=zPi[a]*(*q); q++;
      }
      if(z1==0){error++;}
      else{
	zinew[a]=1./z1;
	float d=zi[a]-zinew[a]; dsumi+=d*d; zi[a]=zinew[a];
      }
    }
    for(a=0; a<Naa; a++){
      if(zjnew[a]==0){error++;}
      else{
	zjnew[a]=1./zjnew[a];
	float d=zj[a]-zjnew[a]; dsumj+=d*d; zj[a]=zjnew[a];
      }
    }
    if(error){
      printf("ERROR in Normalize_Qprime, %d values are zero\n", error);
      break;
    }

    float e=dsumi+dsumj;
    if(e<e_min){
      e_min=e; for(a=0; a<Naa; a++){ziopt[a]=zi[a]; zjopt[a]=zj[a];}
      if(e_min < EPS)break;
    }
  }
  if(e_min > EPS){
    int Na2=Naa*Naa;
    printf("WARNING, the normalization of Q did not converge di=%.2g dj=%.2g\n",
	   dsumi, dsumj);
    printf("i_c= %d j_c= %d cont=%d ij=%.0f\n",
	   pair->i1->i_c, pair->i2->i_c, cont, pair->ij);
    if(0){
      printf("Q_prime= ");
      for(a=0; a<Na2; a++){printf(" %.2g", Qprime[a]);} printf("\n");
    }
    printf("P_i= ");
    for(a=0; a<Naa; a++){printf(" %.2g", pair->P1i1_obs[a]);} printf("\n");
    printf("P_j= ");
    for(a=0; a<Naa; a++){printf(" %.2g", pair->P1i2_obs[a]);} printf("\n");
    if(0){    
      printf("z_i= ");
      for(a=0; a<Naa; a++){printf(" %.2g", zi[a]);} printf("\n");
      printf("z_j= ");
      for(a=0; a<Naa; a++){printf(" %.2g", zj[a]);} printf("\n");
    }
    printf("z_i^opt= ");
    for(a=0; a<Naa; a++){printf(" %.2g", ziopt[a]);} printf("\n");
    printf("z_j^opt= ");
    for(a=0; a<Naa; a++){printf(" %.2g", zjopt[a]);} printf("\n");
    if(e_min>EPS_MAX)return(-1);
  }
  // Final prediction of Q
  float *Q=Q_pred; q=Qprime;
  for(a=0; a<Naa; a++){
    for(b=0; b<Naa; b++){
      *Q=(*q)*(ziopt[a]*zjopt[b]);
      q++; Q++;
    }
  }
  return(e_min);
}

float Normalize_Qprime_symm(float *Q_pred, float *Qprime,
			    struct pair_class  *pair, int cont)
{
  // i1=i2 => zi=zj
  int IT_MAX=500; float EPS=0.0002, EPS_MAX=0.01, e_min=1000;
  int a, b, iter;
  double z1, zij=0, z[Naa], zP[Naa], zopt[Naa];
  // Initialization, 
  for(a=0; a<Naa; a++){z[a]=1;}
  float *P1i1=pair->P1i1_obs, *q=Qprime;
  // Normalize Qprime
  for(a=0; a<Naa; a++){
    z1=0;
    for(b=0; b<Naa; b++){z1+=(*q)*P1i1[b]; q++;}
    zij+=P1i1[a]*z1;
  }
  zij=fabs(zij);
  if(zij<=0 || isnan(zij)){
    printf("ERROR 0 in Normalize_Qp_symm, zij= %.2g\n", zij);
    return(-1);
  }
  zij=1./zij;
  if(zij<=0 || isnan(zij)){
    printf("ERROR 1 in Normalize_Qp_symm, zij= %.2g\n", zij);
    return(-1);
  }

  q=Qprime;
  for(a=0; a<Naa; a++){
    for(b=0; b<Naa; b++){
      if(isnan(*q)){
	printf("ERROR in Normalize_Qp_symm, q(%d,%d) is nan\n",a,b);
	return(-1);
      }
      (*q)*=zij; q++;
    }
  }

  for(iter=0; iter<IT_MAX; iter++){
 
    float e=0; int error=0;
    for(a=0; a<Naa; a++){zP[a]=z[a]*P1i1[a];}
    q=Qprime;
    for(a=0; a<Naa; a++){
      z1=0;
      for(b=0; b<Naa; b++){z1+=(*q)*zP[b]; q++;}
      if(z1==0 || isnan(z1)){error++;}
      else{
	z1=1./z1;
	float d=z[a]-z1; e+=d*d; z[a]=z1;
      }
    }
    if(error){
      printf("ERROR in Normalize_Qprime_symm iter %d "
	     "%d values are zero or nan\n", iter, error);
      return(-1);
    }
    if(e<e_min){
      e_min=e; for(a=0; a<Naa; a++){zopt[a]=z[a];}
    }
    if(e < EPS)break;
  }
  if(e_min>EPS){
    int Na2=Naa*Naa;
    printf("WARNING, the normalization of Q (symmetric) did not converge"
	   " e_min=%.2g\n", e_min);
    printf("i_c= %d j_c= %d cont=%d ij=%.0f\n",
	   pair->i1->i_c, pair->i2->i_c, cont, pair->ij);
    printf("Q_prime= ");
    for(a=0; a<Na2; a++){printf(" %.2g", Qprime[a]);} printf("\n");
    printf("P= ");
    for(a=0; a<Naa; a++){printf(" %.2g", pair->P1i1_obs[a]);} printf("\n");
    printf("z= ");
    for(a=0; a<Naa; a++){printf(" %.2g", z[a]);} printf("\n");
    printf("z^opt= ");
    for(a=0; a<Naa; a++){printf(" %.2g", zopt[a]);} printf("\n");
    if(e_min>EPS_MAX)return(-1);
  }
  // Final prediction of Q
  float *Q=Q_pred; q=Qprime;
  for(a=0; a<Naa; a++){
    for(b=0; b<Naa; b++){
      if(zopt[a] && zopt[b]){
	*Q=(*q)*(zopt[a]*zopt[b]);
      }else{
	*Q=1.0;
      }
      q++; Q++;
    }
  }
  return(e_min);
}

float Update_z_naive(float *zi, float *zj, float *Q, struct pair_class *pair)
{

  /* Q_ij(a,b)= Q'_ij(a,b)xi(a)xj(b)
   sum_b Q'_ij(a,b)xi(a)xj(b)P_j(b)=1 forall i same for j =>
   1/xi(a)=sum_b Q'_ij(a,b)xj(b)Pj(b)
   We fix the scale imposing that sum_a xi(a)/20 = 1

   This normalization is necessary because of the invariance
   xi(a)=k*xi(a), xj(a)=(1/k)*xj(a). We choose the gauge so that <xi>=1
   Normalizing both xi and xj makes dKL and lik much worse
   Not normalizing any of them makes lik and dKL a bit better
   but it makes Lambda less robust.
   The best option seems to be to normalize only xi */

  float dsumi=0, dsumj=0; double z;
  float sumi=0, sumj=0, zinew[Naa], zjnew[Naa], pz[Naa], *p, *q; 
  int a, b;

  // Update zj
  for(a=0; a<Naa; a++)pz[a]=zi[a]*pair->P1i1_obs[a];
  for(b=0; b<Naa; b++){
    z=0; q=Q+b; p=pz; for(a=0; a<Naa; a++){z+=(*q)*(*p); q+=Naa; p++;}
    z=1/z; zjnew[b]=z; sumj+=z*pair->P1i2_obs[b];
  }
  if(sumj==0){printf("ERROR j in normalize\n");}
  else{
    for(a=0; a<Naa; a++){
      //zjnew[a]/=sumj;
      float d=fabs(zj[a]-zjnew[a]); dsumj+=d; zj[a]=zjnew[a];
    }
  }

  // Update and normalize zi
  for(a=0; a<Naa; a++)pz[a]=zj[a]*pair->P1i2_obs[a];
  q=Q;
  for(a=0; a<Naa; a++){
    z=0; p=pz; for(b=0; b<Naa; b++){z+=(*q)*(*p); q++; p++;}
    z=1/z; zinew[a]=z; sumi+=z*pair->P1i1_obs[a]; 
  }
  if(sumi==0){printf("ERROR i in normalize\n");}
  else{
    dsumi=0;
    for(a=0; a<Naa; a++){
      zinew[a]/=sumi;
      float d=fabs(zi[a]-zinew[a]); dsumi+=d; zi[a]=zinew[a];
    }
  }
  return(dsumi+dsumj);
}

float Normalize_Qprime_old(float *Q, struct pair_class  *pair)
{
  // Normalize so that sum_a P_ia (Q_iajb-1)=0
  if(pair->n<=NTHR){
    printf("ERROR, nothing to normalize\n"); exit(8);
  }
  int VBS=0; // Verbose

  // Optimization parameters
  int IT_MAX=400, IT1=IT_MAX/2, iter, nini=0, a, b;
  float EPS=0.0004; // Mean increment = EPS/40= 0.00001
  float error, e_thr=ERR_THR; // Mean error = e_thr/30 = 0.0005
  float e_min=1000, e_min_0, e_ini=e_min;

  // Dynamic variables
  float zi[Naa], zj[Naa], ziopt[Naa], zjopt[Naa], yi[Naa], yj[Naa];
  // Initialization: Normalize Q'
  float *q =Q, *p; double sum=0, z;
  //Check_nan_f(Q, Na2, "Normalize Qprime", -1);
  for(a=0; a<Naa; a++){
    z=0; p=pair->P1i2_obs; for(b=0; b<Naa; b++){z+=(*q)*(*p); q++; p++;}
    yi[a]=z; sum+=pair->P1i1_obs[a]*z;
  }
  for(a=0; a<Naa; a++){
    z=0; q=Q+a; p=pair->P1i2_obs;
    for(b=0; b<Naa; b++){z+=(*q)*(*p); q+=Naa; p++;}
    yj[a]=z;
  }
  q=Q;
  //Check_nan_f(Q, Na2, "Normalize Qprime", -2);
  // Initialize zi=1, zj=1
  if(sum && isinf(sum)==0){
    for(a=0; a<Na2; a++){(*q)/=sum; q++;}
    for(a=0; a<Naa; a++){zi[a]=yi[a]/sum; zj[a]=yj[a]/sum;}
    printf("sum= %.2g\n", sum);
    Check_nan_f(Q, Na2, "Normalize Qprime", -3);
  }else{
    printf("ERROR in Normalize_Qprime, sum is %.2g\n", sum);
    printf("prot: L=%.0f U=%.4f pair: ij= %.1f <Cij>= %.2g C_nat= %d  ",
	   pair->prot_class->L, pair->prot_class->U_ave,
	   pair->ij, pair->Cont_freq, pair->C_nat);
    printf("i= %d j= %d\n", pair->i1->i, pair->i2->i);
    printf("P1i_obs: ");
    for(a=0; a<Naa; a++){printf(" %.2f", pair->P1i1_obs[a]);}
    printf("\n");
    printf("P1j_obs: ");
    for(a=0; a<Naa; a++){printf(" %.2f", pair->P1i2_obs[a]);}
    printf("\n");
    printf("Q:\n");
    for(a=0; a<Naa; a++){
      for(b=0; b<Naa; b++){printf(" %.2f", Q[a*Naa+b]);} printf("\n");
    }
    exit(8);
  }
  // Changed, 23/07/2025
  //Check_nan_f(Q, Na2, "Normalize Qprime", -4);

  float d=0;
  for(iter=0; iter<IT_MAX; iter++){

    if(iter<IT1){
      d=Update_z_Newton(zi, zj, yi, yj, Q, pair);
      Check_nan_f(Q, Na2, "Qprime update_Newton", iter);
    }else if(iter==IT1){
      e_min_0=e_min;
      for(a=0; a<Naa; a++){zi[a]=1; zj[a]=1;}
      d=Update_z_naive(zi, zj, Q, pair);
      Check_nan_f(Q, Na2, "Qprime update_naive", iter);
    }else{
      d=Update_z_naive(zi, zj, Q, pair);
      Check_nan_f(Q, Na2, "Qprime update_naive", iter);
    }

    // Compute error
    float qnew[Na2];
    Update_q(qnew, Q, zi, zj);
    error=Error_norm(qnew, pair);
    if(VBS)printf("e= %.4g\n", error);   
    if(error<e_min){
      e_min=error; for(a=0; a<Naa; a++){ziopt[a]=zi[a]; zjopt[a]=zj[a];}
    }
    if(error<e_thr)break;
   
    // Change zk if no convergence
    if(d<EPS){
      int it=iter; if(it>=IT1)it-=IT1;
      float x=(float)it/IT1, f=20/x; int k=2*Naa*x;
      if(VBS)
	printf("z converges but Q does not. Changing z %d factor %.2g\n", k,f);
      if(nini==0){e_ini=e_min;} nini++;
      for(a=0; a<Naa; a++){zi[a]=1; zj[a]=1;}
      if(k<Naa){zj[k]*=f;}else{zi[k-Naa]*=f;}
    }

  }
  // Use optimal z always
  for(a=0; a<Naa; a++){zi[a]=ziopt[a]; zj[a]=zjopt[a];}
  
  // Final prediction of Q
  // xi = zi/Pi
  Update_q(Q, Q, zi, zj);
  Check_nan_f(Q, Na2, "Q Update_q", -1);

  // Test normalization condition
  error=Error_norm(Q, pair);
  if(error > e_thr){
    printf("ERROR, Q_pred is not well normalized err= %.4g\n", error);
    printf("%d initializations e_ini= %.4g err_Newton= %.4g\n",
	   nini, e_ini, e_min_0);
    printf("prot: L=%.0f U=%.4f pair: ij= %.1f <Cij>= %.2g C_nat= %d  ",
	   pair->prot_class->L, pair->prot_class->U_ave,
	   pair->ij, pair->Cont_freq, pair->C_nat);
    printf("i= %d j= %d\n", pair->i1->i, pair->i2->i);
    float qmax=Q[0], qmin=qmax; int c;
    for(c=1; c<Na2; c++){
      if(Q[c]>qmax){qmax=Q[c];}else if(Q[c]<qmin){qmin=Q[c];}
    }
    printf("Q: Max= %.2g Min= %.2g\n", qmax, qmin);
    if(VBS && (printconv==0)){
      Print_error_Q(pair, Q, zi,zj, ziopt, zjopt);
    }
    printconv++;
  } // end of test 
  return(error);
}




float Update_z_Newton(float *zi, float *zj, float *yi, float *yj,
		      float *Q, struct pair_class *pair)
{
  /* Q_ij(a,b)= Q'_ij(a,b)xi(a)xj(b)
   F_ia == X_ia sum_b Q'_ia,jb P_jb X_jb  - 1  == X_ia*Y_ia -1
   dF_ia/dX_ic = delta_ac Y_ia
   dF_ia/dX_jb = X_ia Q'_ia,jb P_jb
   0=F_ia(X0)~ F_ia(X)+ dF_ia/dX_ia*(X0_ia-X_ia)+sum_b dF_ia/dX_jb*(X0_jb-X_jb)
   = (X_ia*Y_ia-1)+Y_ia*(X0_ia-X_ia)+X_ia*sum_b Q'_ia,jb P_jb*(X0_jb-X_jb)
   = X0_ia*Y_ia -  X_ia*Y_ia + X_ia*Y_ia(X0_j) -1 =>
   X0_ia = X_ia - (X_ia*Y_ia(X0_j)-1)/Y_ia(X_j) 
  */
  Check_nan_f(Q, Na2, "Qprime update_Newton", -1);

  float dsumi=0, dsumj=0; double z; int a, b;
  float sumi=0, sumj=0, zinew[Naa], zjnew[Naa], yinew[Naa], yjnew[Naa];
  float pz[Naa], *p, *q, d; 
  
  // Update zj
  for(a=0; a<Naa; a++)pz[a]=zi[a]*pair->P1i1_obs[a];
  for(b=0; b<Naa; b++){
    z=0; q=Q+b; p=pz; for(a=0; a<Naa; a++){z+=(*q)*(*p); q+=Naa; p++;}
    yjnew[b]=z;
    zjnew[b]=zj[b];
    if(yj[b]){zjnew[b]-=(zj[b]*z-1)/yj[b];}
    sumj+=zjnew[b]*pair->P1i2_obs[b];
  }
  Check_nan_f(Q, Na2, "Qprime update_Newton", -2);
  if(sumj==0){printf("ERROR j in normalize\n");}
  else{
    for(a=0; a<Naa; a++){
      //zjnew[a]/=sumj;
      d=fabs(zj[a]-zjnew[a]); dsumj+=d; zj[a]=zjnew[a]; yj[a]=yjnew[a];
    }
  }
  
  // Update and normalize zi
  for(a=0; a<Naa; a++)pz[a]=zj[a]*pair->P1i2_obs[a];
  q=Q;
  for(a=0; a<Naa; a++){
    z=0; p=pz; for(b=0; b<Naa; b++){z+=(*q)*(*p); q++; p++;}
    yinew[a]=z;
    zinew[a]=zi[a];
    if(yi[a]){zinew[a]-=(zi[a]*yinew[a]-1)/yi[a];}
    sumi+=zinew[a]*pair->P1i1_obs[a];
  }
  Check_nan_f(Q, Na2, "Qprime update_Newton", -3);
  if(sumi==0){printf("ERROR i in normalize\n");}
  else{
    for(a=0; a<Naa; a++){
      //zinew[a]/=sumi;
      d=fabs(zi[a]-zinew[a]); dsumi+=d; zi[a]=zinew[a]; yi[a]=yinew[a];
    }
  }
  return(dsumi+dsumj);
}

void Update_q(float *qnew, float *Q, float *zi, float *zj)
{
  int a, b; float *q=qnew, *qold=Q; 
  for(a=0; a<Naa; a++){
    for(b=0; b<Naa; b++){
      (*q)=(*qold)*zi[a]*zj[b]; q++; qold++;
    }
  }
}

float Error_norm(float *Q, struct pair_class *pair)
{
  float err_i=0, err_j=0, *q; int a, b;
  for(a=0; a<Naa; a++){
    double P1=0; q=Q+a*Naa; float *p=pair->P1i1_obs;
    for(b=0; b<Naa; b++){P1+=(*q)*(*p); q++; p++;}
    P1-=1; err_i+=fabs(P1);
    P1=0; q=Q+a; p=pair->P1i2_obs;
    for(b=0; b<Naa; b++){P1+=(*q)*(*p); q+=Naa; p++;}
    P1-=1; err_j+=fabs(P1);
  }
  return(err_i+err_j);
}

void Print_error_Q(struct pair_class *pair, float *Q,
		   float * zi, float * zj, float * ziopt, float * zjopt)
{
  int a;
  //printf("Phi0= %.2g Phi01= %.2g Phi11= %.2g Phi2= %.2g\n",
  //	 pair->Phi0, pair->Phi01, pair->Phi11, pair->Phi2);
  printf("Q_prime= ");
  for(a=0; a<Na2; a++){printf(" %.2g", Q[a]);} printf("\n");
  printf("P_i= ");
  for(a=0; a<Naa; a++){printf(" %.3f", pair->P1i1_obs[a]);} printf("\n");
  printf("P_j= ");
  for(a=0; a<Naa; a++){printf(" %.3f", pair->P1i2_obs[a]);} printf("\n");
  printf("z_i= ");
  for(a=0; a<Naa; a++){printf(" %.2g", zi[a]);} printf("\n");
  printf("z_j= ");
  for(a=0; a<Naa; a++){printf(" %.2g", zj[a]);} printf("\n");
  printf("z_i^opt= ");
  for(a=0; a<Naa; a++){printf(" %.2g", ziopt[a]);} printf("\n");
  printf("z_j^opt= ");
  for(a=0; a<Naa; a++){printf(" %.2g", zjopt[a]);} printf("\n");
  Check_nan_f(Q, Na2, "Q in print error", 0);
}

/******************** Auxiliary routines *********************/

float **Energy_over_T(float **E_cont_gap, int Naa, float TEMP)
{
  float **E_cont_T=malloc(Naa*sizeof(float *));
  for(int a=0; a<Naa; a++){
    E_cont_T[a]=malloc(Naa*sizeof(float));
    for(int b=0; b<Naa; b++)E_cont_T[a][b]=E_cont_gap[a][b]/TEMP;
  }
  return(E_cont_T);
}

float Average_energy_over_T(int *i_aa, int L, float TEMP)
{
  // Computes the average contact energy for a given sequence

  // Normalize contact frequencies: not necessary

  //float norm[Naa];
  double E_ave=0, Zc=0;
  for(int l=IJ_MIN; l<L; l++){
    float cij=Cont_f[l]; //cij=1;
    for(int i=0; i<(L-l); i++){
      int a=i_aa[i]; if(a<0 && Naa==20){continue;} // gap
      int b=i_aa[i+l]; if(b<0 && Naa==20){continue;}
      E_ave+=cij*E_cont_T[a][b]; Zc+=cij;
    }
  }
  //float sT=sqrt(TEMP);
  return(E_ave/(TEMP*Zc));
}

float Corr_coeff(float *offset, float *slope, float *xx, float *yy, int n){
  double x1=0, x2=0, y1=0, y2=0, xy=0;
  int i; float *x=xx, *y=yy;
  for(i=0; i<n; i++){
    x1 += *x; x2+= (*x)*(*x);
    y1 += *y; y2+= (*y)*(*y);
    xy += (*x)*(*y); x++; y++;
  }
  float nx2=n*x2-x1*x1, ny2=n*y2-y1*y1;
  *slope=(n*xy-y1*x1); if(nx2)(*slope)/=nx2;
  *offset=(y1-(*slope)*x1)/n;
  float r=(*slope);
  if((nx2)&&(ny2))r*=sqrt(nx2/ny2);
  return(r);
}


float Corr_coeff_w(float *offset, float *slope,
		   float *xx, float *yy, double *ww, int n)
{
  double x1=0, x2=0, y1=0, y2=0, xy=0, w1=0;
  int i; float *x=xx, *y=yy;
  double *w=ww;
  for(i=0; i<n; i++){
    double wx=(*w)*(*x), wy=(*w)*(*y);
    x1 += wx; x2+= wx*(*x);
    y1 += wy; y2+= wy*(*y);
    xy += wx*(*y); w1+=*w;
    x++; y++; w++;
  }

  double nx2=w1*x2-x1*x1, ny2=w1*y2-y1*y1;
  *slope=(w1*xy-y1*x1);
  if(nx2>0){(*slope)/=nx2;}else{nx2=0;}
  *offset=(y1-(*slope)*x1)/w1;
  float r=(*slope); 
  if((nx2>=0)&&(ny2>0)){
    r*=sqrt(nx2/ny2);
  }else{
      printf("ERROR in Corr_coeff, negative variance: ");
      printf("X2= %.2g Y2= %.2g W1= %.2g\n", nx2/w1, ny2/w1, w1);
      return(-100);
  }
  return(r);
}

void Copy_vec(float *v1, float *v2, int N){
  for(int a=0; a<N; a++)v1[a]=v2[a];
}

void Fill_C_mat(int **C, int nres, struct contact *cont_list, int Nc){
  int i, j; struct contact *cl=cont_list;
  for(i=0; i<nres; i++)for(j=0; j<nres; j++)C[i][j]=0;
  for(i=0; i<Nc; i++){
    C[cl->res1][cl->res2]=1; C[cl->res2][cl->res1]=1; cl++; 
  }
}

int Find_p_class(int nres, float U_ave,
		 int N_L, int *L_bin, int N_U, float *U_bin)
{
  // Place L and then U
  int i=Find_bin_i(nres, N_L, L_bin);
  int j=Find_bin_f(U_ave,N_U, U_bin);
  return(i*N_U+j);
}

int Find_ij_class(int ij, int C_mat, int C2, int N_ij, int *ij_bin,
		  int i1, int i2, int N_C, int N_C2)
{
  // Place ij and then C_nat
  int i=Find_bin_i(ij, N_ij, ij_bin),j;
  int ncp;  // =(N_C*i1)-i1+(i2-i1)=(N_C-2)*i1+i2 
  if(i1==0){ncp=i2;}
  else if(i1==1){ncp=2+i2;}
  else{ncp=5;}
 

  if(C_mat){j=1;}
  else if(C2){j=2;}
  else{j=0;}
  int n=ncp*3*N_ij+j*N_ij+i; // 6=N_C*(N_C-1)/2
  return(n);
}

struct pair_class *Find_pair(struct pair_class *pair_c, int Npair_prot,
			     int i1, int i2, int C,
			     int ij, int N_ij, int *ij_bin)
{
  // Place ij and then C_nat
  int ij_index=Find_bin_i(ij, N_ij, ij_bin);
  //int C=C_mat; if(C2 && C==0){C=2;} // Contact: 0, 1, 2
  struct pair_class *pair=pair_c;
  for(int j=0; j<Npair_prot; j++){
    if(pair->ij_index==ij_index &&
       pair->C_nat==C &&
       pair->i1->i==i1 &&
       pair->i2->i==i2){
      return(pair);
    }
    pair++;
  }
  printf("ERROR, pair i1= %d i2= %d C= %d ij_index= %d "
	 "could not be found among these %d pairs\n"
	 "# i1 i2 C ij\n", i1, i2, C, ij_index, Npair_prot);
  pair=pair_c;
  for(int j=0; j<Npair_prot; j++){
    printf("%d\t%d\t%d\t%d\n",
	   pair->i1->i, pair->i2->i, pair->C_nat, pair->ij_index);
    pair++;
  }
  exit(8);
  
  return(NULL);
}

int Find_bin_i(int L, int N_L, int *L_bin){
  int i; for(i=0; i<N_L; i++)if(L<=L_bin[i])return(i);
  return(N_L-1);
}

int Find_bin_f(float L, int N_L, float *L_bin){
  int i; for(i=0; i<N_L; i++)if(L<=L_bin[i])return(i);
  return(N_L-1);
}

/*********************** Observed Q *****************************/

double Marginalize_N2_symm(float *P1_i, double *N2_obs)
{
  //float *P1_i=pair->P1i1_obs, *P1_j=pair->P1i2_obs;
  //double *N2_obs=pair->N2_obs; 
  double n=0; int a, b;
  for(a=0; a<Naa; a++){
    double P1=0, *N2=N2_obs+a*Naa;
    for(b=0; b<Naa; b++){P1+=*N2; N2++;}
    P1_i[a]=P1;
    n+=P1;
  }
  if(n==0){
    printf("ERROR in Marginalize_N2, no data found\n"); exit(8);
    return(0);
  }
  float Psumi=0;
  for(a=0; a<Naa; a++){
    P1_i[a]/=n;
    Psumi+=P1_i[a];
  }
  if(fabs(Psumi-1.0)>0.001){
    printf("ERROR, wrongly normalized P1: %.2g n= %.4g\n",
	   Psumi, n);
    for(a=0; a<Naa; a++){printf("%.3g\n", P1_i[a]);}
    exit(8);
  }
  return(n);
}

void Pairwise_Q(float *Q_obs, float *P1i1, float *P1i2, double *N2_obs)
{
  // n=sum_c N2[c]
  int a, b;
  double nn=0;
  float Pi[Naa], Pj[Naa];
  for(a=0; a<Naa; a++){Pi[a]=0; Pj[a]=0;}
  double *N2=N2_obs;
  for(a=0; a<Naa; a++){
     for(b=0; b<Naa; b++){
       Pi[a]+=*N2;
       Pj[b]+=*N2;
       N2++;
     }
     nn+=Pi[a];
  }
  for(a=0; a<Naa; a++){
    Pi[a]/=nn; P1i1[a]=Pi[a];
    Pj[a]/=nn; P1i2[a]=Pj[a];
  }

  float *Q=Q_obs; N2=N2_obs;
  for(a=0; a<Naa; a++){
    float na=Pi[a]*nn;
    for(b=0; b<Naa; b++){
      if((Pi[a]>0)&&(Pj[b]>0)){
	*Q=*N2/(na*Pj[b]);
      }else{
	if(*N2!=0.0){
	  printf("ERROR: Wrong probabilities P1i= %.4f P1j= %.4f",
		 Pi[a], Pj[b]);
	  printf(" N2= %.0f\n", *N2); exit(8);
	}
	*Q=1.0;
      }
      Q++; N2++;
    }
  }
  Check_nan_f(Q_obs, Na2, "Q_obs in Pairwise_Q", 0);
}

/************************** Main computation **************************/

float Score(float Lambda, struct prot_class *prot, int Npair_prot, 
	    int I_SCORE, int comp_all)
{
  double Score=0, norm=0;
  int indirect=0;
  // indirect contacts are computed only if comp_all

  Check_nan_f(prot->Q_global, Na2, "Q_global in Score", -1);

  /****************** Predict Q for contacts, no contacts, indirect *********/

  int N_ind=0;
  for(int j=0; j<Npair_prot; j++){
    struct pair_class *pair=prot->pairs+j;
    float w;
    if(comp_all){w=pair->n;}else{w=pair->w;}
    if(pair->n<=NTHR)continue;
    if(pair->C_nat==2 && comp_all==0)continue;

    pair->C_pred=0;
    for(int cont=0; cont<=2; cont++){

      float *Q_pred=pair->Q_pred[cont],
	*log_Q_pred=pair->log_Q_pred[cont];
      if(cont<=1){
	// Q'(a,b)=Phi1*U(a,b)+Phi2*U(a,b)^2
	Compute_Q(Q_pred, Lambda, pair, cont);
	Compute_logarithm(log_Q_pred, Q_pred, qmin, lmin, Na2);
      }else{ // cont =2
	if(comp_all==0){break;}
	if(PDB){
	  N_ind+=Indirect_log_Q_PDB(log_Q_pred, pair, prot->pairs, Npair_prot);
	  for(int c=0; c<Na2; c++){Q_pred[c]=exp(log_Q_pred[c]);}
	}else{
	  Indirect_Q_ali(Q_pred, j, prot->pairs);
	  Compute_logarithm(log_Q_pred, Q_pred, qmin, lmin, Na2);
	}
      }
 
      /****************** Compute score, sum it over pairs ****************/

      char what[90];
      sprintf(what, "log(Q) Lambda=%.2f cont=%d pair c1=%d c2=%d ij=%.1f",
	      Lambda, cont, pair->i1->i_c, pair->i2->i_c, pair->ij);
      Check_nan_f(log_Q_pred, Na2, what, cont);

      pair->score.log_lik[cont]=Likelihood(log_Q_pred, pair);  
      if(isnan(pair->score.log_lik[cont])){
	printf("ERROR, log.lik is nan cont=%d Lambda= %.2g\n", cont, Lambda);
      }      
      
      // Kullback-Leibler divergence
      if(comp_all || I_SCORE==1){
	pair->score.dKL[cont]=dKL(log_Q_pred, pair);
      }
      // Relative Mean Error
      if(comp_all || I_SCORE==2){
	pair->score.MSE[cont]=MSE(log_Q_pred, pair);
      }
      
      if((I_SCORE==0 &&
	  pair->score.log_lik[cont] > pair->score.log_lik[pair->C_pred]) ||
	 (I_SCORE==1 && 
	  pair->score.dKL[cont] < pair->score.dKL[pair->C_pred]) ||
	 (I_SCORE==2 && 
	  pair->score.MSE[cont] < pair->score.MSE[pair->C_pred]))
	{
	  pair->C_pred=cont;
	}
    } // end cont

    // Score to be optimized: exclude indirect contacts if indirect==0
    if(pair->C_nat==2 && indirect==0){continue;}
    if(I_SCORE==0){
      Score+= w*pair->score.log_lik[pair->C_nat];
    }else if(I_SCORE==1){
      Score-= w*pair->score.dKL[pair->C_nat];
    }else{
      Score-= w*pair->score.MSE[pair->C_nat];
    }
    norm+=w;  
  } // end pairs
  if(comp_all){
    printf("Indirect contacts computed for %d over %d pairs\n",
	   N_ind, 3*Npair_prot);
  }

  if(norm==0){
    printf("ERROR, no data found for computing the score\n");
    exit(8);
  } 
  return(Score/norm);
}

float **Compute_E_cont_gap(int Naa){
  if((Naa!=20)&&(Naa!=21)){
    printf("ERROR, Naa=%d but instead of 20 a.a. or 21 (20+gap)\n", Naa);
    exit(8);
  }
  // Opening a gap is scored as the minimum energy of an a.a.
  float **E=Allocate_mat2_f(Naa, Naa), E_gap[Naa]; int a,b;
  for(a=0; a<20; a++){
    E_gap[a]=0;
    for(b=0; b<20; b++){
      E[a][b]=Econt[a][b];
      if(E[a][b]<E_gap[a])E_gap[a]=E[a][b];
    }
  }
  if(Naa==21){
    float E_gg=100;
    for(a=0; a<20; a++){
      E[20][a]=-E_gap[a]; E[a][20]=-E_gap[a];
      if((a==0)||(E_gap[a]>E_gg))E_gg=E_gap[a];
    }
    E[20][20]=-E_gg;
  }
  return(E);
}

/************************* Indirect contacts *******************/

void Indirect_Q_ali(float *Q, int i_ij, struct pair_class *pairs)
{
  if(label_ij==NULL){
    printf("ERROR, no label available for pairs of aligned residues\n");
    exit(8);
  }
  struct pair_class *pair_ij=pairs+i_ij;
  int i=pair_ij->i1->i, j=pair_ij->i2->i, c;
  double Ind[Na2], Qk[Na2]; for(c=0; c<Na2; c++)Ind[c]=0;
  for(int k=0; k<pair_ij->prot_class->L; k++){
    if((k==i)||(k==j))continue;
    int ik, jk;
    if(k>i){ik=label_ij[i][k];}else{ik=label_ij[k][i];}
    if(k>j){jk=label_ij[j][k];}else{jk=label_ij[k][j];}
    if((ik<0)||(jk<0))continue;
    Indirect_Q_k(Qk, Naa, pair_ij, pairs+ik, pairs+jk);
    for(c=0; c<Na2; c++)Ind[c]+=Qk[c];
  }
  float Qprime[Na2];
  for(c=0; c<Na2; c++){
    if(isnan(Ind[c])){
      printf("ERROR in Indirect_Qali, c=%d\n", c); exit(8);
    }
    Qprime[c]=Ind[c];
  }
  Normalize_Qprime(Q, Qprime, pair_ij, 2);
}

int Indirect_log_Q_PDB(float *log_Q, struct pair_class *pair_ij,
			struct pair_class *pairs, int Npair_prot)
{
  // Find most likely intermediate contact
  int n_ind=0;
  float d_thr=4, log_lik_opt=-10000;
  float Q_ind[Na2], Q_prime[Na2], log_Q_ind[Na2];

  for(int ik=0; ik<Npair_prot; ik++){

    struct pair_class *pair_ik=pairs+ik;
    if(pair_ik->C_nat!=1 || pair_ik->n < 10 )continue;

    if(pair_ij->i1 != pair_ik->i1 &&
       pair_ij->i1 != pair_ik->i2 && 
       pair_ij->i2 != pair_ik->i1 &&
       pair_ij->i2 != pair_ik->i2){continue;}

    for(int jk=ik; jk<Npair_prot; jk++){

      struct pair_class *pair_jk=pairs+jk;
      if(pair_jk->C_nat!=1 || pair_jk->n < 10 )continue;
      
      if(Indirect_Q_new(Q_prime, pair_ij, pair_ik, pair_jk, d_thr)<0)continue;

      if(1){
	float *Qij0=pair_ij->Q_pred[0];
	for(int c=0; c<Na2; c++){Q_prime[c]*=Qij0[c];}
      }

      if(Normalize_Qprime(Q_ind, Q_prime, pair_ij, 2)<0){continue;}
      Compute_logarithm(log_Q_ind, Q_ind, qmin, lmin, Na2);

      float log_lik=Likelihood(log_Q_ind, pair_ij);
      if(n_ind==0 || log_lik>log_lik_opt){
	log_lik_opt=log_lik;
	Copy_vec(log_Q, log_Q_ind, Na2);
      }
      n_ind++;
    }
  }

  if(n_ind==0){
    printf("WARNING, no candidate indirect contacts found. "
	   "n1=%d n2=%d ij=%.1f, set no contact\n",
	   pair_ij->i1->i_c, pair_ij->i2->i_c, pair_ij->ij);
    Copy_vec(log_Q, pair_ij->log_Q_pred[0], Na2);
  }
  return(n_ind);
}

int Match_site(struct pair_class *pair_ij, struct pair_class *pair_ik,
	       int *i1, int *i2)
{
  int nk=0;
  if(pair_ij->i1==pair_ik->i1){
    i1[nk]=1; i2[nk]=-1; nk++;
  }else if(pair_ij->i1==pair_ik->i2){
    i1[nk]=2; i2[nk]=-1; nk++;
  }
  if(pair_ij->i1==pair_ij->i2){return(nk);}
  if(pair_ij->i2==pair_ik->i1){
    i1[nk]=-1; i2[nk]=1; nk++;
  }else if(pair_ij->i2==pair_ik->i2){
    i1[nk]=-1; i2[nk]=2; nk++;
  }
  return(nk);
}

int Indirect_Q_new(float *Q, struct pair_class *pair_ij,
		   struct pair_class *pair_1, struct pair_class *pair_2,
		   float d_thr)
{
  // sites 1,2: position 1 and 2 of pair_1
  // sites 3,4: position 1 and 2 of pair_2
  struct site *i1=pair_ij->i1, *i2=pair_ij->i2;
  struct site *s1=pair_1->i1, *s2=pair_1->i2, *s3=pair_2->i1, *s4=pair_2->i2;

  // Check sites
  // Q[a][b]=sum_c Qa[ca]*Qb[cb]*Pk[c]
  // ca=fa*a*fca*c;
  // cb=fb*b+fcb*c;
  // There are only 6 possibilities, since s1>=s2 and s1>=s3>=s4
  // and i1 and i2 cannot be in the same pair
  float *Q1, *Q2, *Pk1, *Pk2;
  int f1=0, f2=0, fc1=0, fc2=0; // c
  if(i1==s1 && i2==s3){
    if(s2!=s4){return(-1);}
    Q1=pair_1->Q_pred[1]; Pk1=pair_1->P1i2_obs;
    Q2=pair_2->Q_pred[1]; Pk2=pair_2->P1i2_obs;
    f1=Naa; fc1=1; // Qa[c1]=Q1[a][c] c1=a*Naa+c;
    f2=Naa; fc2=1; // Q2[c2]=Q2[b][c] c2=b*Naa+c;
  }else if(i1==s1 && i2==s4){
    if(s2!=s3){return(-1);}
    Q1=pair_1->Q_pred[1]; Pk1=pair_1->P1i2_obs;
    Q2=pair_2->Q_pred[1]; Pk2=pair_2->P1i1_obs;
    f1=Naa; fc1=1; // Q1[c1]=Q1[a][c] c1=a*Naa+c;
    f2=1; fc2=Naa; // Qb[c2]=Q2[c][b] c2=c*Naa+b;
  }else if(i1==s2 && i2==s3){
    if(s1!=s4){return(-1);}
    Q1=pair_1->Q_pred[1]; Pk1=pair_1->P1i1_obs;
    Q2=pair_2->Q_pred[1]; Pk2=pair_2->P1i2_obs;
    f1=1; fc1=Naa; // Q1[c1]=Q1[c][a] c1=c*Naa+a;
    f2=Naa; fc2=1; // Q2[c2]=Q2[b][c] c2=b*Naa+c;
  }else if(i1==s2 && i2==s4){
    if(s1!=s3){return(-1);}
    Q1=pair_1->Q_pred[1]; Pk1=pair_1->P1i1_obs;
    Q2=pair_2->Q_pred[1]; Pk2=pair_2->P1i1_obs;
    f1=1; fc1=Naa; // Q1[c1]=Q1[c][a] c1=c*Naa+a;
    f2=1; fc2=Naa; // Q2[c2]=Q2[c][b] c2=c*Naa+b;
  }else if(i1==s3 && i2==s2){
    if(s1!=s4){return(-1);}
    Q1=pair_2->Q_pred[1]; Pk1=pair_2->P1i2_obs;
    Q2=pair_1->Q_pred[1]; Pk2=pair_1->P1i1_obs;
    f1=Naa; fc1=1; // Q1[c1]=Q1[a][c] c1=a*Naa+c;
    f2=1; fc2=Naa; // Q2[c2]=Q2[c][b] c2=c*Naa+b;
  }else if(i1==s4 && i2==s2){
    if(s1!=s3){return(-1);}
    Q1=pair_2->Q_pred[1]; Pk1=pair_2->P1i1_obs;
    Q2=pair_1->Q_pred[1]; Pk2=pair_1->P1i1_obs;
    f1=1; fc1=Naa; // Q1[c1]=Q1[c][a] c1=c*Naa+a;
    f2=1; fc2=Naa; // Q2[c2]=Q2[c][b] c2=c*Naa+b;
  }else{
    return(-1);
  }
  
  // Check range
  if(pair_ij->ij >= pair_1->ij && pair_ij->ij >= pair_2->ij){
    // k between i and j, ij= ik+ jk
    if(fabs(pair_ij->ij-(pair_1->ij+pair_2->ij))>d_thr){return(-1);}
  }else{
    // k outside ij, ij=|ik-jk|
    float d=fabs(pair_2->ij-pair_1->ij)-pair_ij->ij;
    if(fabs(d)>d_thr){return(-1);} //d<0 ||
  }
  
  float *q=Q;
  for(int a=0; a<Naa; a++){
    int a1=f1*a;
    for(int b=0; b<Naa; b++){
      double QQ=0; int c1=a1, c2=f2*b;
      for(int c=0; c<Naa; c++){
	QQ+=(Pk1[c]+Pk2[c])*(Q1[c1]-1)*(Q2[c2]-1);
	c1+=fc1; c2+=fc2;
      }
      *q=1+0.5*QQ; q++;
    }
  }
  
  return(0);
}

float Compute_indirect(float *Q, struct pair_class *pair, double *P)
{
  float Qprime[Na2];
  for(int c=0; c<Na2; c++){
    if(isnan(P[c])){
      printf("ERROR in Compute_indirect! c=%d\n", c); exit(8);
    }
    Qprime[c]=P[c]-1;
  }
  float e=Normalize_Qprime(Q, Qprime, pair, 2);
  return(e);
}

void Indirect_Q_k(double *P_ind, int Naa,
		  struct pair_class *pair_ij,
		  struct pair_class *pair_ik,
		  struct pair_class *pair_jk)
{
  double *q=P_ind;
  float *Qik=pair_ik->Q_pred[1], *Qjk=pair_jk->Q_pred[1];
  int a, b, c;

  if((pair_ij->ij>=pair_ik->ij)&&(pair_ij->ij>=pair_jk->ij)){
    // i < k < j
    float *P1k=pair_ik->P1i2_obs;
    for(a=0; a<Naa; a++){
      int na=a*Naa;
      for(b=0; b<Naa; b++){
	double Q=0; int ac=na, cb=b;  
	for(c=0; c<Naa; c++){
	  Q+=P1k[c]*Qik[ac]*Qjk[cb];
	  ac++; cb+=Naa;
	}
	(*q)=Q; q++;
      }
    }
  }else if((pair_ik->ij<pair_jk->ij)){ // k < i < j
    float *P1k=pair_ik->P1i1_obs;
    for(a=0; a<Naa; a++){
      for(b=0; b<Naa; b++){
	double Q=0; int ca=a, cb=b;	  
	for(c=0; c<Naa; c++){
	  Q+=P1k[c]*Qik[ca]*Qjk[cb];
	  ca+=Naa; cb+=Naa;
	}
	(*q)=Q; q++;
      }
    }
  }else{ // i < j < k
    float *P1k=pair_ik->P1i2_obs;
    for(a=0; a<Naa; a++){
      int na=Naa*a;
      for(b=0; b<Naa; b++){
	double Q=0; int ac=na, bc=b*Naa;	  
	for(c=0; c<Naa; c++){
	  Q+=P1k[c]*Qik[ac]+Qjk[bc]; ac++; bc++;
	}
	(*q)=Q; q++;
      }
    }
  }

}

void Indirect_contacts(int **C2, int **k_indirect, int **l_indirect,
		       struct contact **indirect_cont, int *Nc_indirect,
		       int nres, int **C1)
{
  /* Returns the matrix C2 of contacts through an indirect residue that
     are not direct contacts */
  int ic=0;
  int i, j, k;
  for(i=0; i<nres; i++){
    int *Ci=C1[i];
    for(j=i+1; j<nres; j++){
      int Csum=0, *Cj=C1[j], lmax=0, kmax=-1;
      for(k=0; k<nres; k++){
	if((Ci[k])&&(Cj[k])){
	  Csum++;
	  int l1=abs(i-k), l2=abs(j-k);
	  if(l1<l2){if(l1>lmax){lmax=l1; kmax=k;}}
	  else{if(l2>lmax){lmax=l2; kmax=k;}}
	}
      }
      C2[i][j]=Csum; C2[j][i]=Csum;
      if(Csum){
	k_indirect[i][j]=kmax; l_indirect[i][j]=lmax;
	k_indirect[j][i]=kmax; l_indirect[j][i]=lmax;
	if(Ci[j]==0)ic++;
      }
    }
  }
  int Nic=ic; ic=0;
  (*indirect_cont)=malloc(Nic*sizeof(struct contact));
  for(i=0; i<nres; i++){
    for(j=i+1; j<nres; j++){
      if((C2[i][j])&&(C1[i][j]==0)){
	if(ic>=Nic){
	  printf("ERROR, too many indirect contacts\n");
	  printf("ic= %d i= %d j= %d nres= %d\n", ic, i, j, nres);
	  exit(8);
	}
	(*indirect_cont)[ic].res1=i;
	(*indirect_cont)[ic].res2=j; ic++;
      }
    }
  } 
  *Nc_indirect=ic;
}

struct contact *Indirect_contact_list(int *Nc_indirect,
				      struct contact *cont_list,
				      int Nc, int nres, int ij_min)
{
  /* Return the list of contacts that are mediated through one or more
     intermediate residues but are not direct contacts */
  int i, j, *indirect[nres], *direct[nres];
  for(i=0; i<nres; i++){
    indirect[i]=malloc(nres*sizeof(int));
    direct[i]=malloc(nres*sizeof(int));
    for(j=0; j<nres; j++){
      indirect[i][j]=0; direct[i][j]=0;
    }
  }
  struct contact *c1=cont_list;
  for(int ic=0; ic<Nc; ic++){
    direct[c1->res1][c1->res2]=1;
    direct[c1->res2][c1->res1]=1;
    c1++;
  }
  c1=cont_list;
  for(int ic=0; ic<Nc; ic++){
    struct contact *c2=c1+1;
    for(int jc=ic+1; jc<Nc; jc++){
      if(c1->res1==c2->res1){
	i=c1->res2; j=c2->res2;
      }else if(c1->res2==c2->res1){
	i=c1->res1; j=c2->res2;
      }else if(c1->res2==c2->res2){
	i=c1->res1; j=c2->res1;
      }else{
	goto next;
      }
      if(abs(i-j)<ij_min)goto next;

      indirect[i][j]=1;
      indirect[j][i]=1;
    next:
      c2++;
    }
    c1++;
  }

  *Nc_indirect=0;
  for(i=0; i<nres; i++){
    for(j=i+ij_min; j<nres; j++){
      if(indirect[i][j]==0)continue;
      if(direct[i][j]){indirect[i][j]=0;}
      else{(*Nc_indirect)++;}
    }
  }
  struct contact *ind_cont=malloc((*Nc_indirect)*sizeof(struct contact));
  struct contact *ind=ind_cont; int ic=0;
  for(i=0; i<nres; i++){
    for(j=i+ij_min; j<nres; j++){
      if(indirect[i][j]){
	ind->res1=i;
	ind->res2=j;
	ind++;
	ic++;
      }
    }
  }
  for(i=0; i<nres; i++){
    free(indirect[i]);
    free(direct[i]);
  }
  return(ind_cont);
}

/************************ Scores *******************************/

void Compute_logarithm(float *log_Q, float *Q, float qmin, float lmin, int Na2)
{
  float *q=Q, *lq=log_Q;
  for(int c=0; c<Na2; c++){
    if(*q>qmin){*lq=log(*q);}else{*lq=lmin;} q++; lq++;
  }
}

float MSE(float *log_Q_pred, struct pair_class *pair)
{
  double d2=0, dd=0; //d1=0, norm=0; 
  for(int c=0; c<Na2; c++){
    float w=pair->N2_obs[c];
    if(w==0)continue;
    dd+=w*pow((pair->log_Q_obs[c]-log_Q_pred[c]),2);
    d2+=w*pow(pair->log_Q_obs[c], 2);
  }
  if(d2){dd/=d2;}
  return(dd);
}

float Likelihood(float *log_Q_pred, struct pair_class *pair)
{
  double logP=0;
  for(int c=0; c<Na2; c++){
    logP+=pair->N2_obs[c]*log_Q_pred[c];
  }
  return(logP/pair->n);
}

float dKL(float *log_Q_pred, struct pair_class *pair)
{
  double d=0;
  for(int c=0; c<Na2; c++){
    d+=pair->N2_obs[c]*(pair->log_Q_obs[c]-log_Q_pred[c]);
  }
  return(d/pair->n);
}

float Compute_specific_mut_inf(struct pair_class *pair){
  double M2=0;
  float *log_Q_glob=pair->prot_class->log_Q_glob;
  for(int c=0; c<Na2; c++){
    if(pair->N2_obs[c]){
      M2+=pair->N2_obs[c]*(pair->log_Q_obs[c]-log_Q_glob[c]); 
    }
  }
  return(M2/pair->n);
}

float Global_mut_inf(struct pair_class *pair, int Na2){
  double M2=0;
  float *log_Q_glob=pair->prot_class->log_Q_glob;
  for(int c=0; c<Na2; c++){
    if(pair->N2_obs[c])M2+=pair->N2_obs[c]*log_Q_glob[c];
  }
  return(M2/pair->n);
}

void Contact_probabilities(struct pair_class *pairs, int Npair_prot)
{
  /* float C_EXP=0.5; // P(C_ij=1)=<C_ij>^C_EXP, 0<=C_EXP<=1
  // int APC=1; // Perform Average Pair Correction on d_lik?
  // float SD=S_DEPEND; // Number of dependent sequences */

  // Likelihood difference for all pairs and Cij
  int i; struct pair_class *pair=pairs;
  double P_cont_pred=0, P_cont_obs=0, f1;
  for(i=0; i<Npair_prot; i++){
    pair=pairs+i;
    pair->score.d_log_lik=(pair->score.log_lik[1]-pair->score.log_lik[0]);
    if(pair->score.d_log_lik < 0){
      float f=exp(pair->score.d_log_lik);
      f1=f/(1+f);
    }else{
      f1=1./(1+exp(-pair->score.d_log_lik));
    }
    float cij=pair->Cont_freq, fc=f1*cij;;
    pair->score.P_C_nat=fc/(fc+(1-cij));
    //pair->score.P_C_nat=cij/(cij+ff*(1-cij));
    P_cont_pred+=pair->n*pair->score.P_C_nat;
    P_cont_obs+=pair->n*pair->Cont_freq;
  }
  printf("Computing contact probabilities.\n");

  // Normalize so that predicted and observed contacts are equal?
  if(0){
    double normP=P_cont_obs/P_cont_pred;
    if((P_cont_pred > 0) && (normP>0)){
      printf("Correction factor: %.2g\n",normP);
      for(i=0; i<Npair_prot; i++)pairs[i].score.P_C_nat*=normP;
    }
  }

  // Important: reduce prob if E_cont >0 or Mut_inf<0
  double Mut_inf_1=0, Mut_inf_2=0;
  double E_1=0, E_2=0;
  for(i=0; i<Npair_prot; i++){
    pair=pairs+i;
    Mut_inf_1+=pair->mut_inf;
    Mut_inf_2+=pair->mut_inf*pair->mut_inf;
    E_1+=pair->E_cont_obs;
    E_2+=pair->E_cont_obs*pair->E_cont_obs;
  }
  E_1/= Npair_prot; E_2-=Npair_prot*E_1*E_1;
  E_2=sqrt(E_2/(Npair_prot-1));
  Mut_inf_1/= Npair_prot; Mut_inf_2-= Npair_prot*Mut_inf_1*Mut_inf_1;
  Mut_inf_2=sqrt(Mut_inf_2/(Npair_prot-1));
  float E_thr=E_1+E_2, M_thr=Mut_inf_1+Mut_inf_2;

  pair=pairs; 
  for(i=0; i<Npair_prot; i++){
    if(0){
      if(pair->E_cont_obs > E_thr)
	pair->score.P_C_nat *= exp(-(pair->E_cont_obs-E_thr)/E_2);
      if(pair->mut_inf<M_thr)
	pair->score.P_C_nat *= exp((pair->mut_inf-M_thr)/Mut_inf_2);
    }
    if(pair->C_pred!=2){
      //if(pair->score.P_C_nat>0.5){pair->C_pred=1;}
      //else{pair->C_pred=0;}
      if(pair->score.d_log_lik > 0){pair->C_pred=1;}
      else{pair->C_pred=0;}
    }
    pair++;
  }

}

void Prot_results(struct prot_class *prot_class, int Npair_prot, int L_tar)
{
  //prot_class->Lambda=Lambda;
  //prot_class->f_indir=f_indir;
  //int Na2=Naa*Naa;

  float E_cont_not[Na2], *E=E_cont_not;
  for(int a=0; a<Naa; a++)for(int b=0; b<Naa; b++){*E=E_cont_gap[a][b]; E++;}
  for(int j=0; j<Npair_prot; j++){
    struct pair_class *pair=prot_class->pairs+j;
    if(pair->n<=NTHR){continue;}
    pair->score.r_op=
      Corr_coeff_w(&(pair->score.offset_op), &(pair->score.slope_op),
		   pair->log_Q_pred[pair->C_pred], pair->log_Q_obs,
		   pair->N2_obs, Na2);

    Zero_mean(E_cont_not, E_cont_1, pair->P1i1_obs, pair->P1i2_obs, Naa);
    pair->cont_score=Contact_score(pair, E_cont_not, Na2);
    pair->score.r_oe=
      Corr_coeff_w(&(pair->score.offset_oe), &(pair->score.slope_oe),
		   pair->log_Q_obs, E_cont_not,
		   pair->N2_obs, Na2);

  }
  // L_tar is used only here
  //APC_pairs(prot_class->pairs, Npair_prot, L_tar, "mut_inf");
  if(prot_class->pairs->C_nat<0)return;
  for(int c=0; c<3; c++){
    prot_class->r_op[c]=0;
    prot_class->r_oe[c]=0;
    prot_class->cont_score[c]=0;
    prot_class->lik[c]=0;
    prot_class->dKL[c]=0;
    prot_class->MSE[c]=0;
    prot_class->norm[c]=0;
  }
  for(int j=0; j<Npair_prot; j++){
    struct pair_class *pair=prot_class->pairs+j;
    if(pair->n<=NTHR){continue;}
    float w=pair->n; // pair->w;
    prot_class->r_op[pair->C_nat]+= w*pair->score.r_op;
    prot_class->r_oe[pair->C_nat]+= w*pair->score.r_oe;
    prot_class->cont_score[pair->C_nat]+= w*pair->cont_score;
    prot_class->lik[pair->C_nat]+= w*pair->score.log_lik[pair->C_nat];
    prot_class->dKL[pair->C_nat]+= w*pair->score.dKL[pair->C_nat];
    prot_class->MSE[pair->C_nat]+= w*pair->score.MSE[pair->C_nat];
    prot_class->norm[pair->C_nat]+=w;
  }
  for(int c=0; c<3; c++){
    if(prot_class->norm[c]==0)continue;
    prot_class->r_op[c]/=prot_class->norm[c];
    prot_class->r_oe[c]/=prot_class->norm[c];
    prot_class->dKL[c]/=prot_class->norm[c];
    prot_class->lik[c]/=prot_class->norm[c];
    prot_class->MSE[c]/=prot_class->norm[c];
    prot_class->cont_score[c]/=prot_class->norm[c];
  }
}

float Contact_score(struct pair_class *pair, float *E_cont_not, int Na2)
{
  double score=0;
  for(int c=0; c<Na2; c++){
    score+=pair->N2_obs[c]*pair->log_Q_obs[c]*E_cont_not[c];
  }
  if(pair->n)score/=pair->n;
  return(-score);
}

/******************* Pairwise statistics ***************************/

void Pairwise_statistics(struct prot_class *prot, int Npairs)
{
  lmin=log(qmin);

  // Observed Q
  int empty[3], c; for(c=0; c<3; c++){empty[c]=0;} 
  printf("Computing observed and specific mutual information\n");
  for(int i=0; i<Npairs; i++){
    struct pair_class *pair=prot->pairs+i;
    if(pair->n<=NTHR){empty[pair->C_nat]++; continue;} 
    Pairwise_Q(pair->Q_obs, pair->P1i1_obs, pair->P1i2_obs, pair->N2_obs);
    // Q_ij(a,b)=P_ij(a,b)/P_i(a)P_j(b)
    Compute_logarithm(pair->log_Q_obs, pair->Q_obs, qmin, lmin, Na2);
  }
  int empty_all=0; int np=Npairs/3;
  for(c=0; c<3; c++){
    empty_all+=empty[c];
    prot->Npair[c]=np-empty[c];
  } 
  //prot->Npairs=(Npairs-empty_all);
  printf("%d classes with data and %d empty\n",Npairs-empty_all, empty_all);
  if(empty_all==Npairs){printf("No classes left, exiting\n"); return;}
}



void Compute_Q_global(struct prot_class *prots, int Np_class,
		      int Npairs)
{
  // global Q computed with optimal regularization
  printf("Computing global propensities Q_global for %d prot classes\n",
	 Np_class);
  double N2_global[Na2]; int a;
  for(a=0; a<Na2; a++)N2_global[a]=0;

  for(int ip=0; ip<Np_class; ip++){
    struct pair_class *pair=prots[ip].pairs;
    for(int i=0; i<Npairs; i++){
	// Exclude short range loops: local optimum at 33 only slightly better 
	// if(pair->ij < 33)continue; 
      for(a=0; a<Na2; a++)N2_global[a]+=pair->N2_obs[a];
      pair++;
    }
  }
  // Symmetrize N2_global
  Symmetrize_mat(N2_global);

  float *P1=prots->P1_glob;
  //double n=
  //Marginalize_N2_symm(P1, N2_global);
  Pairwise_Q(prots->Q_global, P1, P1, N2_global);
  // Q_global(a,b)= N2(a,b)/P1(a)P1(b)

  // Regularize the pairwise probability
  // L1-L2 regularization with maximum specific heat
  double one_third=1.0/(3.0), two_third=2.0*one_third;
  for(a=0; a<Naa; a++){P1[a]=two_third*P1[a]+one_third;}
  if(1){
    for(a=0; a<Na2; a++){
      prots->Q_global[a]=two_third*prots->Q_global[a]+one_third;
    }
  }else{
    for(a=0; a<Na2; a++){
      prots->Q_global[a]=two_third*(prots->Q_global[a]-1);
    }
    Zero_mean(prots->Q_global, prots->Q_global, P1, P1, Naa);
    for(a=0; a<Na2; a++){prots->Q_global[a]+=1;}
  }

  Check_nan_f(prots->Q_global, Na2, "Q_global in Compute Q_global", 0);
  Compute_logarithm(prots->log_Q_glob, prots->Q_global, qmin, lmin, Na2);
  Copy_vec(prots->Q_glob_ini, prots->Q_global, Na2);

  for(int ip=1; ip<Np_class; ip++){
    Copy_vec((prots+ip)->Q_glob_ini, prots->Q_global, Na2);
    Copy_vec((prots+ip)->Q_global, prots->Q_global, Na2);
    Copy_vec((prots+ip)->log_Q_glob, prots->log_Q_glob, Na2);
    Copy_vec((prots+ip)->P1_glob, prots->P1_glob, Naa);
  }
}

void Symmetrize_mat(double *X)
{
  // It is not correct to symmetrize, unless the site is the same
  for(int a=0; a<Naa; a++){
    int ab=a*Naa, ba=a;
    for(int b=0; b<a; b++){
      float xx=0.5*(X[ab]+X[ba]);
      X[ab]=xx; X[ba]=xx; ab++; ba+=Naa;
    }
  }
}

int Update_Q_global(struct prot_class *prots, int Np_class, int Npair_prot)
{
  // global Q
  printf("Updating unspecific mutual information Q_global\n");

  double Q_pred_sum[Na2], norm=0; int a;
  for(a=0; a<Na2; a++){Q_pred_sum[a]=0;}

  for(int ip=0; ip<Np_class; ip++){
    struct pair_class *pairs=prots[ip].pairs;
    for(int i=0; i<Npair_prot; i++){
      struct pair_class *pair=pairs+i;
      if(pair->C_nat==2){continue;}
      if(pair->n<=NTHR)continue;
      float w=pair->n; // pair->w
      if(isnan(w) || isinf(w)){continue;}
      float *QQ=pair->Q_pred[pair->C_nat];
      for(a=0; a<Na2; a++){Q_pred_sum[a]+=w*QQ[a];}
      norm+=w;
    }
  }
  if(norm==0){
    printf("ERROR in Update_Q_global, w= %.1g\n", norm); return(-1);
  } 
  int nan=Check_nan(Q_pred_sum, Na2, "Q_glob_update", -1);
  if(nan){return(-1);}

  for(a=0; a<Na2; a++){
    if(Q_pred_sum[a]){
      prots->Q_global[a]=(prots->Q_glob_ini[a]*prots->Q_global[a])
	/(Q_pred_sum[a]/norm);
    }else{
      prots->Q_global[a]=1;
    }
  }
  //Check_nan(prots->Q_global, Na2, "Qnew_global", 3);
  Compute_logarithm(prots->log_Q_glob, prots->Q_global, qmin, lmin, Na2);

  for(int ip=1; ip<Np_class; ip++){
    Copy_vec((prots+ip)->Q_global, prots->Q_global, Na2);
    Copy_vec((prots+ip)->log_Q_glob, prots->log_Q_glob, Na2);
  }
  return(0);

}

int Check_nan(double *Q, int n, char *what, int step)
{
  double *q=Q; int error_nan=0, error_inf=0;
  for(int a=0; a<n; a++){
    if(isnan(*q)){error_nan++;}
    else if(isinf(*q)){error_inf++;}
  }

  if(error_nan || error_inf){
    printf("ERROR, %s is nan %d times and inf %d times.",
	   what, error_nan, error_inf);
    if(step>=0){printf(" step=%d", step);}
    printf("\n");
    return(-1);
    //exit(8);
  }
  return(0);
}

int Check_nan_f(float *Q, int n, char *what, int step)
{
  for(int a=0; a<n; a++){
    if(isnan(Q[a])){
      printf("ERROR, %s is nan step=%d a=%d\n", what, step, a);
      exit(8);
      //return(-1);
    }
    if(isinf(Q[a])){
      printf("ERROR, %s is inf step=%d a=%d\n", what, step, a);
      exit(8);
      //return(-1);
    }
  }
  return(0);
}


void APC_pairs(struct pair_class *pairs, int Npairs, int L, char *what)
{
  // Compute Average Product Correction (APC) to mutual information or P_cont
  printf("Computing Average Product Correction to %s\n", what);
  int i_what, i;
  if(strncmp(what, "mut_inf", 7)==0){i_what=0;}
  else if(strncmp(what, "d_lik", 5)==0){i_what=1;}
  else{printf("ERROR in APC, undefined %s\n", what); exit(8);}

  double average[L], norm[L];
  for(i=0; i<L; i++){
    average[i]=0;
    norm[i]=0;
  } 
  struct pair_class *pair=pairs; float wy;
  for(i=0; i<Npairs; i++){
    if(i_what==0){
      wy=pair->mut_inf;
    }else{
      wy=pair->w*exp(pair->score.log_lik[1]-pair->score.log_lik[0]);
    }
    if(pair->i1->i_c<L){
      average[pair->i1->i_c]+=wy;
      norm[pair->i1->i_c]+=pair->w;
    }
    if(pair->i2->i_c<L){
      average[pair->i2->i_c]+=wy;
      norm[pair->i2->i_c]+=pair->w;
    }
    pair++;
  }

  double ave_all=0, norm_all=0;
  for(i=0; i<L; i++){
    if(norm[i]){
      ave_all+=average[i];
      norm_all+=norm[i];
      average[i]/=norm[i];
      if(i_what)average[i]=log(average[i]);
    }
  }
  if(norm_all)ave_all/=norm_all;
  if(i_what)ave_all=log(ave_all);
  for(i=0; i<Npairs; i++){
    pair=pairs+i;
    float APC=0; int i1=pair->i1->i, i2=pair->i2->i;
    if(i1>=0 && i1<L && i2>=0 && i2<L){
      APC=average[i1]*average[i2];
    }
    if(ave_all)APC/=ave_all;
    if(i_what==0){
      pair->mut_inf-=APC;
    }else{
      pair->score.d_log_lik-=APC;
    }
  }
}

void Print_Cont_stat(int L, int ij_min)
{
  char name_out[100]; sprintf(name_out, "Cont_stat_%d.dat", L);
  FILE *file_out=fopen(name_out, "w");
  printf("Writing contact statistics in file %s\n", name_out);

  float Cont_norm[L], nc_nc_norm[L], CNc1_norm[L], CNc2_norm[L], Cnc1_norm[L];
  Normalize_cont_freq(Cont_norm, Cont_f,
		      Cnc1_norm, Cnc1,
		      nc_nc_norm, nc_nc,
		      CNc1_norm, Cz1,
		      CNc2_norm, Cz2,
		      NULL, NULL,
		      Nc1L[L], Nc2L[L], Nc3L[L], 0, L);

  fprintf(file_out, "#|i-j| 2=<Cij> 3=<Cij*Nc>-<Cij><Nc> 4=<ni*nj>-<ni*nj>");
  fprintf(file_out, " 5=<Cij*Nc^2>-<Cij><Nc^2>-2(<Cij*Nc>-<Cij><Nc>)");
  fprintf(file_out, " 6=<Cij*ni>-<Cij><ni>\n");
  for(int l=ij_min; l<L; l++){
    fprintf(file_out, "%d\t%.4f", l, Cont_norm[l]);
    fprintf(file_out, "\t%.4g", CNc1_norm[l]);
    fprintf(file_out, "\t%.4g", nc_nc_norm[l]);
    fprintf(file_out, "\t%.4g", CNc2_norm[l]);
    fprintf(file_out, "\t%.4g", Cnc1_norm[l]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}


