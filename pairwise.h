struct score{
  float P_C_nat; // Probability of contact
  float log_lik[3]; // no cont, cont, indirect cont
  float d_log_lik;
  float dKL[3];
  float MSE[3];
  float r_op, slope_op, offset_op; // log(Q) observed vs predicted  
  float r_oe, slope_oe, offset_oe; // log(Q) observed vs energy
};

// Program: CNc1, CNc1_norm

struct site{
  struct prot_class *prot_class;
  int i; // index of site, for grouped sites it is equal to i_c
  int n; // number of instances
  float nc; // Number of contacts (mean)
  int i_c; // index of number of contacts
  float P1_obs[21]; // Frequency  of amino acid at the site
  float U1[21];     // U1_i[a]=sum_b P1_i[b]U[a][b]
  float cU1[21];    // cU1_i[a]=sum_j <C_ij>U1_j[b]
};

struct pair_class{
  // Class description:
  int C_nat; // Native contacts 0,1,2
  int ij_index; // Contact range
  //int i_c, j_c; // index of number of contacts, now in site
  // Pointers
  struct prot_class *prot_class;
  struct site *i1, *i2;
  // Measured quantities
  double n; // Number of sampled pairs
  double ij;   // Average contact range
  float w; // weight: 0, n, or sqrt/n)
  float E_cont_obs; // Mean energy of sampled pairs
  float mut_inf;   // Mutual information
  float min_dist; // Average of minimal distance between pairs 
  double *N2_obs; // Observed pairwise amino acid distribution
  float P1i1_obs[21], P1i2_obs[21];
  float *Q_obs;
  float *log_Q_obs;
  //float nc1, nc2; // Number of contacts of residues 1 and 2, now in site
  // Measured quantities for making predictions
  //float Phi1, Phi11, Phi2; // outdated
  double Cont_freq;     // <C_ij>
  double Cont_freq_Nc;  // <C_ij*Nc>-<C_ij><Nc>
  double Cont_freq_Nc2; // <C_ij*Nc^2>-<C_ij><Nc^2>
  double Cont_freq_nc;  // <C_ij*n_i>-<C_ij><n_i>
  double nc_nc;         // <ni*nj>-<ni>*<nj>
  double num_cont;
  double Cont_freq_indir;
  // Predictions
  double cont_score;
  int C_pred; // Most likely type of contact 0,1,2
  float **ene; //E[0][ab]: no contact E[1][ab]: contact
  //float **ene1, **ene2, *ene12;
  //float bij; // <Cij>/(1-<Cij>)
  //float logb;
  float *log_Q_pred[3]; // Prediciton for no contact, contacts or indirect cont
  float *Q_pred[3]; 
  //float *Q_pred_cont; 
  //float *Q_pred_pred; 
  //float *P2_null;
  // Scores
  struct score score;
  // indirect contacts
  int ik, jk;
  float P_ind;
  /* long lk; 
  float *P1k;
  int ik_indirect;
  int jk_indirect;
  int k_internal; */
  float error;
};

struct prot_class{
  // Class description:
  int i_L;
  double L; // length
  double Nc; // Number of contacts
  int i_U;  // Average energy
  double U_ave;
  float Lambda; // selection parameter
  float P1_glob[21];
  float *Q_global;
  float *log_Q_glob;
  float *Q_glob_ini;
  float *Q_glob_opt;
  //float U1[21];
  double U_norm;
  //float U1ave_overT[441]; // 21*21
  int Np; // Number of proteins
  float f_indir;
  //float mu; // Normalization of log-likelihood
  int npairs; // Number of pairs
  int Npair[3];
  double cont_score[3];
  double r_op[3];
  double r_oe[3];
  double lik[3];
  double dKL[3];
  double MSE[3];
  double norm[3];
  //float *Phi1;
  //float *Phi2;
  //float *Cont_freq;
  //float *Cont_freq_indir;
  //float *Phi3;
  struct site *site;
  struct pair_class *pairs;
};


