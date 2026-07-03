/*****************************
 *         Structures        *
 *****************************/
struct contact{
  short res1, res2;
};
struct protein{
  char name[10];
  char chain[10];
  char exp_meth;
  int oligomer;
  int nres;
  int *i_aa;
  char *amm;
  int Nc;
  short **C_mat;
  struct contact *cont_list;
  struct contact *indirect_contacts;
  float **min_dist;
  int Nc_indirect;
  int *n_cont;
  float U_ave;
  float U1[20]; 
};







