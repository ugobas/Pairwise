#include "coord.h"
#include "protein4.h"
#include "read_pdb.h"
//#include "codes.h"
#include "allocate.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
extern int Verbose;

#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWH-"
#define ASTPATH    "/data/ortizg/databases/astral_40/"
#define PDBPATH    "/data/ortizg/databases/pdb/"
#define PDBEXT     ".pdb" /*.ent.Z*/
#define PDBCAT     "zcat"  /* zcat */
#define PDBTMP     "pdb.dat"

extern char dir_pdb[200];
static char res_exo[400][5], res_std[400][5];
static int n_atom;
#define ATOM_MAX 7000000
#define RES_MAX 10000
atom atom_read[ATOM_MAX];
struct residue res[RES_MAX];
int init_map=0;
float cont_thr_c=4.5, cont_thr_c2;

static int Next_residue(int *N_res, int *start, struct residue *seq,
			atom *first_atom, short *i_atom,
			char *res_type_old, int *res_num_old, char *icode_old,
			char *chain_old, int *hetatm_old, int n_exo,
			char *res_type, int res_num, char icode, char chain,
			int hetatm);
static short Write_residue(char *res_type_old, atom *first_atom, short i_atom,
			   struct residue *ptr_tmp, int n_exo, int res_num,
			   char icode, char chain, int hetatm);
static char Code_3_1(char *res);
int Code_AA(char res);
int Get_compression(char *pdb_name);
void Get_name(char *name, char *file_name);
static struct contact *
Compute_cont_list_2(int *n_cont, int *Nc, float **min_dist,
		    struct residue *res, int N_res, int ij_min);
static int Contact(struct residue *res_i, struct residue *res_j);
static int Count_proteins_seq(FILE *file_seq);
float **Compute_min_dist(struct residue *res, int N_res);
static float Min_dist(struct residue *res_i, struct residue *res_j);
int Get_NMR(char *pdb_file);
int Get_monomeric(char *pdb_file);
void Set_conts_new( struct contact **cont, int *nc,
		    int res_old, int res_last, int res_max,
		    int ij_min, int ij_min_new,
		    struct protein *prot);

int Read_proteins_PDB(struct protein **prot_ptr, int N_prot,
		      char **file_pdb, char **chain,
		      char *dir_pdb, char *ext, int ij_min,
		      int MONOMERIC, int Xray)
{
  int NP_MAX=N_prot*20;
  *prot_ptr=malloc(NP_MAX*sizeof(struct protein *));
  struct protein *prot=*prot_ptr;
  struct residue res[L_MAX];
  int n_pdb=0, i_pdb, i;
  int ANISOU, nmr;
  
  for(i_pdb=0; i_pdb<N_prot; i_pdb++){
    char file[400];
    sprintf(file, "%s%s%s", dir_pdb, file_pdb[i_pdb], ext);
    int N_mod=Count_models_PDB(file), imod, kmod;
    if(N_mod<0){
      printf("WARNING, no PDB file found: %s\n", file);
      continue;
    }
    printf("PDB %s, %d models\n", file_pdb[i_pdb], N_mod);
    if(MONOMERIC && (Get_monomeric(file)==0))continue;
    if(Xray && Get_NMR(file))continue;

    for(imod=1; imod<=N_mod; imod++){
      if(N_mod==1){kmod=-1;}else{kmod=imod;}
      strcpy(prot->chain, chain[i_pdb]);
      prot->nres=
	Read_coord(file, &nmr, res, atom_read, prot->chain, &ANISOU, kmod);

      if(prot->nres<=0){
	printf("WARNING, no chain found: %s chain %s\n",
	       file, prot->chain); continue;
      }
      Get_name(prot->name, file_pdb[i_pdb]);
      strcat(prot->name, prot->chain);
      prot->amm=malloc(prot->nres*sizeof(char));
      prot->i_aa=malloc(prot->nres*sizeof(short));
      prot->n_cont=malloc(prot->nres*sizeof(int));
      for(i=0; i<prot->nres; i++){
	if((res[i].i_aa<0)||(res[i].i_aa>=20)){
	  printf("Error, wrong amino acid code %d\n", res[i].i_aa);
	  exit(8);
	}
	prot->i_aa[i]=res[i].i_aa;
	prot->amm[i]=AMIN_CODE[prot->i_aa[i]];
      }
      prot->indirect_contacts=NULL;
      /*prot->cont_list=
	Compute_cont_list(prot->n_cont, &prot->Nc, res, prot->nres, ij_min);
	prot->min_dist=Compute_min_dist(res, prot->nres);*/
      prot->min_dist=Allocate_mat2_f(prot->nres, prot->nres);
      prot->cont_list=
	Compute_cont_list_2(prot->n_cont, &prot->Nc, prot->min_dist,
			    res, prot->nres, ij_min);

      prot++; n_pdb++;
      if(n_pdb>=NP_MAX){
	printf("ERROR, too many input proteins (more than %d)\n", NP_MAX);
	exit(8);
      }
    }
  }

  return(n_pdb);
}

int Read_processed_proteins(struct protein **prot_ptr,
			    char *FILE_STR, char *FILE_SEQ,
			    int IJ_MIN, int IJ_MIN_OLD,
			    int MONOMERIC, int Xray, int L_min)
{
  int IJ_ADD=IJ_MIN_OLD-IJ_MIN;
  FILE *file_str=fopen(FILE_STR, "r");
  if(file_str==NULL){
    printf("ERROR, file %s with contact matrices not found\n",
	   FILE_STR); exit(8);
  }
  FILE *file_seq=fopen(FILE_SEQ, "r");
  if(file_seq==NULL){
    printf("ERROR, file %s with sequences not found\n",
	   FILE_SEQ); exit(8);
  }

  // Count proteins and allocate
  int N_prot=Count_proteins_seq(file_seq);
  *prot_ptr=malloc(N_prot*sizeof(struct protein));

  // Read sequences
  fclose(file_seq);
  file_seq=fopen(FILE_SEQ, "r");
  printf("Reading %s\n", FILE_SEQ);
  char string[1000];
  int n=-1, i; struct protein *prot=*prot_ptr-1;
  fgets(string, sizeof(string), file_seq);
  while(fgets(string, sizeof(string), file_seq)!=0){
    if(string[0]=='#'){
      prot++; n++; i=0;
      if(n > N_prot){
	printf("ERROR file %s more than %d proteins found\n",
	       FILE_SEQ, N_prot); exit(8);
      }
      char experiment[10]; int oligo;
      int m=sscanf(string+1, "%s%d%s%d",
		   prot->name, &prot->nres, experiment, &oligo);
      if(m>2){prot->exp_meth=experiment[0];}
      else{prot->exp_meth=' ';}
      if(m>3){prot->oligomer=oligo;}
      else{prot->oligomer=-1;}
      prot->amm=malloc(prot->nres*sizeof(char));
      prot->i_aa=malloc(prot->nres*sizeof(int));
      prot->n_cont=malloc(prot->nres*sizeof(int));
      for(int k=0; k<prot->nres; k++){prot->n_cont[k]=0;}
    }else{
      char aa=string[0]; if(aa=='*'){aa='A';}
      int a=Code_AA(aa);
      if(a==-2){
	printf("WARNING no AA read in %s: i=%d L=%d %s\n",
	       string, i, prot->nres, prot->name);
      }
      //if(a<0){if(Naa==21){a=20;}else{a=0;}}
      prot->i_aa[i]=a;
      prot->amm[i]=aa;
      i++;
      if(i>prot->nres){
	printf("ERROR protein %s file %s more than %d residues found\n",
	       prot->name, FILE_SEQ, prot->nres); exit(8);
      }
    }
  }
  fclose(file_seq);

  // Read contact matrices
  int nres, res1, res2, nc=0; char name[10];
  struct contact *cont;
  fgets(string, sizeof(string), file_str);
  n=-1; prot=*prot_ptr-1;
  int Nc, res_old=-1;
  while(fgets(string, sizeof(string), file_str)!=0){
    if(string[0]=='#'){
      if(IJ_ADD && n>=0){
	Set_conts_new(&cont,&nc,res_old,nres-1,nres,IJ_MIN_OLD,IJ_MIN,prot);
	if(nc!= prot->Nc){
	  printf("WARNING %d, few contacts %d %d L=%d\n",n,nc,prot->Nc,nres);
	  prot->Nc=nc;
	} 
      }
      prot++; n++; nc=0; res_old=-1;
      if(n > N_prot){
	printf("ERROR file %s more than %d proteins found\n",
	       FILE_STR, N_prot); exit(8);
      }
      sscanf(string+1, "%d%d%s", &nres, &Nc, name);
      prot->Nc=Nc+IJ_ADD*(nres-1)-(IJ_ADD); // Add contacts i,i+1 and i,i+2

      if((nres!=prot->nres)||(strcmp(name, prot->name)!=0)){
	printf("ERROR, I expected protein %s %d ", prot->name, prot->nres);
	printf("but I found protein %s %d\n", name, nres);
	printf("File: %s protein %d\n", FILE_STR, n+1);
	exit(8);
      }
      prot->indirect_contacts=NULL;
      prot->cont_list=malloc(prot->Nc*sizeof(struct contact)); nc=0;
      cont=prot->cont_list;
    }else{
      sscanf(string, "%d%d", &res1, &res2);
      if((res1 > prot->nres)||(res2 > prot->nres)){
	printf("ERROR, read contact %d %d but only %d residues in prot %s\n",
	       res1, res2, prot->nres, prot->name);
	exit(8);
      }
      // New: Set contacts with r+1, r+2
      if(IJ_ADD){
	Set_conts_new(&cont,&nc,res_old,res1,nres,IJ_MIN_OLD,IJ_MIN,prot);
      }
      if(res_old<res1){res_old=res1;}
      if(res2<(res1+IJ_MIN_OLD)){continue;}
      cont->res1=res1; cont->res2=res2;
      prot->n_cont[res1]++;
      prot->n_cont[res2]++;
      nc++; cont++;
      if(nc>prot->Nc){
	printf("ERROR protein %s file %s more than %d contacts found\n",
	       prot->name, FILE_STR, prot->Nc); exit(8);
      }
    }
  }
    if(IJ_ADD){
      Set_conts_new(&cont, &nc, res_old, nres-1, nres,IJ_MIN_OLD,IJ_MIN,prot);
      if(nc!= prot->Nc){
	printf("WARNING %d, few contacts %d %d L=%d\n", n,nc,prot->Nc,nres);
	prot->Nc=nc;
      }
  } 
  fclose(file_str);

  // Check if monomeric and X ray
  n=0;
  int n_short=0, n_NMR=0, n_pol=0;
  printf("Checking monomeric and NMR proteins in %s\n", dir_pdb);
  for(i=0; i<N_prot; i++){
    prot=*prot_ptr+i;
    char file[400]="", pdb[20], *s1=prot->name, *s2=pdb;
    while(*s1!='\0'){if(*s1=='_'){break;} *s2=*s1; s1++; s2++;} *s2='\0';
    if(Xray){
      if(prot->exp_meth!=' '){
	if(prot->exp_meth!='X'){n_NMR++; continue;}
      }else{
	sprintf(file, "%s%s.pdb", dir_pdb, pdb);
	int nmr=Get_NMR(file);
	if(nmr<0){printf("name: %s\n", prot->name); continue;}
	else if(nmr){n_NMR++; continue;}
      }
    }
    if(MONOMERIC){
      if(prot->oligomer>=0){
	if(prot->oligomer!=1){n_pol++; continue;}
      }else{
	sprintf(file, "%s%s.pdb", dir_pdb, pdb);
	int mono=Get_monomeric(file);
	if(mono<0){printf("name: %s\n", prot->name); continue;}
	else if(mono==0){n_pol++; continue;}
      }
    }
    if(prot->nres<L_min){n_short++; continue;}
    if(n!=i)(*prot_ptr)[n]=*prot;
    n++;
  }
  printf("Keeping %d proteins over %d\n", n, N_prot);
  if(n_NMR)printf("%d proteins discarded because NMR\n", n_NMR);
  if(n_pol)printf("%d proteins discarded because not monomeric\n", n_pol);
  if(n_short)printf("%d proteins discarded < %d res\n", n_short, L_min);
  printf("\n");
  N_prot=n;

  return(N_prot);
}

void Set_conts_new( struct contact **cont, int *nc,
		    int res_old, int res_last, int res_max,
		    int ij_min_old, int ij_min,
		    struct protein *prot)
{
  for(int res=res_old+1; res<=res_last; res++){
    for(int r2=res+ij_min; r2<res+ij_min_old; r2++){
      if(r2>=res_max)break;
      (*cont)->res1=res; 
      (*cont)->res2=r2;
      prot->n_cont[res]++;
      prot->n_cont[r2]++;
      (*cont)++; (*nc)++;
    }
  }
  if(*nc>prot->Nc){
    printf("ERROR protein %s (processed) more than %d contacts found\n",
	   prot->name, prot->Nc); exit(8);
  }
}

void Get_name(char *name, char *file_name){
  char *n=name, *s=file_name, *s1=s;
  while(*s!='\0'){if(*s=='/')s1=s; s++;} if(*s1=='/')s1++;
  while(*s1!='\0'){if(*s1=='.'){break;} *n=*s1; n++; s1++;}
  *n='\0';
}

int Read_coord(char *pdb_name, int *nmr, struct residue *seq, atom *atoms,
	       char *chain_to_read, int *ANISOU, int kmod)
{
  int N_res=0, n_exo=0, Compression=0, read_atom=0, start=0;
  FILE *file_in;
  char string[200], command[200];
  atom *atom_ptr=atoms, *first_atom=NULL;
  int N_model=0;

  short i_atom=0, alternative=0;
  int hetatm=0, hetatm_old=0;
  int i, j, res_num, res_num_old=10000;
  char altloc, altloc_sel=' ';
  char chain_old='#', chain='Z';
  char res_type[5], res_type_old[5], icode=0, icode_old;
  float x, y, z;
  char file_name[500];
  int nchains=1, ic, ichain=0, ini_chain=0;
  n_atom=0;

  /* Open file */
  if(Verbose)printf("Reading %s ", pdb_name);
  Compression=Get_compression(pdb_name); 
  if(Compression){
    sprintf(command, "%s %s > %s\n", PDBCAT, pdb_name, PDBTMP);
    system(command); strcpy(file_name, PDBTMP);
  }else{
    sprintf(file_name, "%s", pdb_name);
  }
  file_in=fopen(file_name, "r");
  if(file_in==NULL){
    printf("\nWARNING, file %s not found\n", file_name); return(0);
  }
  //printf("Reading %s\n", file_name);

  // Count chains to read
  for(i=0; i<sizeof(chain_to_read); i++){
    nchains=i;
    if((chain_to_read[i]=='\0')||(chain_to_read[i]==' '))break;
    printf("%c", chain_to_read[i]);    
  }
  if(nchains==0)nchains=1;
  printf(" %d chains to read\n", nchains);

  *nmr=0;
  strcpy(res_type_old,"xxx");
  while(fgets(string, sizeof(string), file_in)!=NULL){

    if(strncmp(string,"ATOM", 4)==0){
      /* Standard residue or DNA basis */
      if((kmod>=0)&&(N_model!=kmod))continue;
      hetatm=0; if(ini_chain==0)ini_chain=1;

    }else if(strncmp(string,"HETATM", 6)==0){
      if((kmod>=0)&&(N_model!=kmod))continue;
      if(ini_chain==0){continue;} hetatm=1;
      /* Cofactor or exotic residue */

    }else if(strncmp(string,"EXPDTA", 6)==0){
      if(strncmp(string+10, "NMR", 3)==0){*nmr=1;} continue;
      /* NMR structure */
    }else if(strncmp(string, "MODEL", 5)==0){
      N_model++; continue;

    }else if((strncmp(string,"TER",3)==0)&&(N_res>0)){
      Next_residue(&N_res, &start, seq, first_atom, &i_atom,
		   res_type_old, &res_num_old, &icode_old,
		   &chain_old, &hetatm_old, n_exo,
		   res_type, res_num, icode, chain, hetatm);
      ini_chain=0;
      continue;

    }else if(strncmp(string,"MODRES", 6)==0){
      res_exo[n_exo][0]=string[12]; res_std[n_exo][0]=string[24];
      res_exo[n_exo][1]=string[13]; res_std[n_exo][1]=string[25];
      res_exo[n_exo][2]=string[14]; res_std[n_exo][2]=string[26];
      for(j=0; j<n_exo; j++){
	if(strncmp(res_exo[j],res_exo[n_exo],3)==0)break;
      }
      if(j==n_exo){n_exo++;} continue;
 
    }else if(strncmp(string,"ENDMDL", 6)==0){
      break;                                    /* end model */
    }else if((strncmp(string,"ANISOU", 6)==0)&&(read_atom)){
      // Anisotropic structure factor
      int i, j;
      if(*ANISOU==0)*ANISOU=1;
      sscanf(string+28, "%f %f %f %f %f %f",
	     &(atom_ptr->anisou[0][0]), &(atom_ptr->anisou[1][1]),
	     &(atom_ptr->anisou[2][2]), &(atom_ptr->anisou[0][1]),
	     &(atom_ptr->anisou[0][2]), &(atom_ptr->anisou[1][2]));
      for(i=0; i<3; i++){
	for(j=i; j<3; j++){
	  atom_ptr->anisou[i][j]*=0.0001;
	  if(j!=i)atom_ptr->anisou[j][i]=atom_ptr->anisou[i][j];
	}
      }
      // Check
      //for(i=0; i<3; i++)B+=aniso[i][i]; B*=26.319; // 8pi^2/2
      //printf("B= %.2f %.2f %d\n", atom_ptr->B_factor, B, n_atom);
      continue;

      /*}else if(strncmp(string,"HELIX ", 6)==0){
	Read_sec_str(string, chain, 'H'); continue;
	}else if(strncmp(string,"SHEET ", 6)==0){
	Read_sec_str(string, chain, 'E'); continue;
	}else if(strncmp(string,"TURN ", 5)==0){
	Read_sec_str(string, chain, 'T'); continue;*/
    }else{
      continue;
    }

    read_atom=0;
    chain=string[21];
    if(*chain_to_read!='*'){
      if((*chain_to_read==' ')||(*chain_to_read=='\0'))*chain_to_read=chain;
      for(ic=0; ic<nchains; ic++)if(chain==chain_to_read[ic])goto read;
      continue;
    }

  read:
    /* Read atom name */
    if((string[13]=='H')||(string[12]=='H') ||
       (string[13]=='D') ||(string[12]=='D'))continue;

    /* Read residue; check if water molecule */
    res_type[0]=string[17]; res_type[1]=string[18]; res_type[2]=string[19];
    res_type[3]='\0';
    if((hetatm==1)&&((strncmp(res_type,"HOH",3)==0)||
		     (strncmp(res_type,"DOD",3)==0)))continue;
    
    icode=string[26]; string[26]=' ';


    /* Read coordinates */
    sscanf(string+22,"%d %f %f %f", &res_num, &x, &y, &z);
    
    /* Check if alternative conformation */
    if((string[72]=='A')&&(string[73]=='L')&&(string[74]=='T')&&
       (string[75]!='1')&&(string[75]!=' '))continue;

    altloc=string[16];
    if(altloc!=' '){
      if(altloc_sel==' '){altloc_sel=altloc;}
      if(altloc!=altloc_sel){continue;}
    }

    if((icode!=icode_old)&&(res_num==res_num_old)){
      if(alternative==1){
	continue;
      }else{
	atom *atom_old=first_atom;
	float dx, dy, dz;
	dx=x-atom_old->r[0]; dy=y-atom_old->r[1]; dz=z-atom_old->r[2];
	if((dx*dx+dy*dy+dz*dz)<.5){
	  alternative=1; continue;  
	}
      }         
    }

    /* New residue */
    if((res_num!=res_num_old)||(icode!=icode_old)||
       (strncmp(res_type, res_type_old, 3)!=0)){
      Next_residue(&N_res, &start, seq, first_atom, &i_atom,
		   res_type_old, &res_num_old, &icode_old,
		   &chain_old, &hetatm_old, n_exo,
		   res_type, res_num, icode, chain, hetatm);
    }
    read_atom=1;
    atom_ptr=atoms+n_atom;
    
    if(i_atom==0)first_atom=atom_ptr;
    i_atom++; n_atom++;
    atom_ptr->r[0] = x;
    atom_ptr->r[1] = y;
    atom_ptr->r[2] = z;
    if(string[12]!=' '){ // Hydrogen atoms
      for(ic=0; ic<4; ic++)atom_ptr->name[ic]=string[12+ic];
    }else{
      for(ic=0; ic<3; ic++)atom_ptr->name[ic]=string[13+ic];
      atom_ptr->name[3]=' ';
    }
    sscanf(string+56, "%f", &atom_ptr->occupancy);
    sscanf(string+60, "%f", &atom_ptr->B_factor);
    atom_ptr->chain=ichain;
  }
  Next_residue(&N_res, &start, seq, first_atom, &i_atom,
	       res_type_old, &res_num_old, &icode_old,
	       &chain_old, &hetatm_old, n_exo,
	       res_type, res_num, icode, chain, hetatm);
  if(N_res >= L_MAX){
    printf("\n ERROR, more than %d residues found\n", L_MAX); exit(8);
  }
  fclose(file_in);
  if(Verbose)printf("%3d residues\n", N_res);
  if(Compression){
    sprintf(command, "rm -f %s\n", PDBTMP); system(command);
  }
  return(N_res);
}

static short Write_residue(char *res_type_old, atom *first_atom, short i_atom,
			   struct residue *ptr_tmp, int n_exo, int res_num,
			   char icode, char chain, int hetatm)
{
  short i, het=1, exo;
  char amm, pdbres[6];

  /* Check amino acid type */
  amm=Code_3_1(res_type_old);
  if(amm != 'X'){
    // Standard residue
    if(hetatm){het=1; goto discard;} // Standard res. and HETATM => cofactor
    het=0; exo=0;
  }else{
    // Non-standard residue
    exo=1;
    for(i=n_exo-1; i>=0; i--){
      if(strncmp(res_type_old,res_exo[i],3)==0){
	amm=Code_3_1(res_std[i]); het=0; break;
      }
    }
  }

  /* Check backbone */
  if(het && (i_atom >=3)){
    // If backbone atoms exist: Modified residue
    int i_N=0, i_CA=0, i_C=0, i;
    atom *atom_ptr=first_atom;
    for(i=0; i<i_atom; i++){
      if(strncmp(atom_ptr->name, "N ", 2)==0){
	i_N=1;
      }else if(strncmp(atom_ptr->name, "CA", 2)==0){
	i_CA=1;
      }else if(strncmp(atom_ptr->name, "C ", 2)==0){
	i_C=1;
      }
      if(i_N && i_CA && i_C){het=0; break;}
      atom_ptr++;
    }
  }

  if(het==0){
    ptr_tmp->atom_ptr=first_atom;
    ptr_tmp->n_atom=i_atom;
    ptr_tmp->amm=amm; ptr_tmp->exo=exo;
    ptr_tmp->i_aa=Code_AA(amm);
    if(ptr_tmp->i_aa<0){
      printf("Unknown residue %s %d%c\n", res_type_old, res_num, icode);
      ptr_tmp->i_aa=0;
    }
    ptr_tmp->chain=chain;
    sprintf(pdbres, "%4d%c", res_num, icode);
    strcpy(ptr_tmp->pdbres, pdbres);
    return(het);
  }

  //printf("%s %d%c  %d %d\n", res_type_old, res_num, icode, hetatm, het);

 discard:

  if(het)
    printf("Group %s %d%c  %c (%d atoms) not a residue\n",
	   res_type_old, res_num, icode, amm, i_atom);

  return(het);
}

int Next_residue(int *N_res, int *start, struct residue *seq,
		 atom *first_atom, short *i_atom,
		 char *res_type_old, int *res_num_old, char *icode_old,
		 char *chain_old, int *hetatm_old, int n_exo,
		 char *res_type, int res_num, char icode, char chain,
		 int hetatm)
{
  if((*start)==0){
    *start=1;
  }else if(*i_atom){
    int het=Write_residue(res_type_old, first_atom, *i_atom, seq+*N_res, n_exo,
			  *res_num_old, *icode_old, *chain_old, *hetatm_old);
    if(het==0){(*N_res)++;}else{n_atom-=(*i_atom);} (*i_atom)=0;
  }
  strcpy(res_type_old,res_type); *icode_old=icode;
  *res_num_old=res_num; *chain_old=chain; *hetatm_old=hetatm;
  
  return(0);
}

char Code_3_1(char *res){
  char code;

  if(strncmp(res,"ALA",3)==0){ code='A';
  }else if(strncmp(res,"GLU",3)==0){ code='E';
  }else if(strncmp(res,"GLN",3)==0){ code='Q';
  }else if(strncmp(res,"ASP",3)==0){ code='D';
  }else if(strncmp(res,"ASN",3)==0){ code='N';
  }else if(strncmp(res,"LEU",3)==0){ code='L';
  }else if(strncmp(res,"GLY",3)==0){ code='G';
  }else if(strncmp(res,"LYS",3)==0){ code='K';
  }else if(strncmp(res,"SER",3)==0){ code='S';
  }else if(strncmp(res,"VAL",3)==0){ code='V';
  }else if(strncmp(res,"ARG",3)==0){ code='R';
  }else if(strncmp(res,"THR",3)==0){ code='T';
  }else if(strncmp(res,"PRO",3)==0){ code='P';
  }else if(strncmp(res,"ILE",3)==0){ code='I';
  }else if(strncmp(res,"MET",3)==0){ code='M';
  }else if(strncmp(res,"PHE",3)==0){ code='F';
  }else if(strncmp(res,"TYR",3)==0){ code='Y';
  }else if(strncmp(res,"CYS",3)==0){ code='C';
  }else if(strncmp(res,"TRP",3)==0){ code='W';
  }else if(strncmp(res,"HIS",3)==0){ code='H';
  }else if(strncmp(res,"HIE",3)==0){ code='H';
  }else if(strncmp(res,"HID",3)==0){ code='H';
  }else if(strncmp(res,"HIP",3)==0){ code='H';
  }else if(strncmp(res,"ASX",3)==0){ code='N';
  }else if(strncmp(res,"GLX",3)==0){ code='Q';
  }else{
    printf("WARNING, a.a. %s not known\n", res);
    return('X');
  }
  return(code);
}

int Code_AA(char res){
  short i; char r;
  i=(int)res; if(i>96){r=(char)(i-32);}else{r=res;}
  for(i=0; i<20; i++)if(r==AMIN_CODE[i])return(i);
  if(res=='X')return(0);
  if((res=='-')||(res=='.')){return(-1);}
  printf("Warning, wrong aa type (%c)\n", res);
  return(-2);
}

int Get_compression(char *pdb_name){
  char *tmp=pdb_name;
  while(*tmp!='\0'){
    if((*tmp=='.')&&(*(tmp+1)=='g')&&(*(tmp+2)=='z'))return(1);
    tmp++;
  }
  return(0);
}

int Count_models_PDB(char *pdb_name){

  int Compression=Get_compression(pdb_name);
  char string[200], command[200], file_name[500];
  if(Compression){
    sprintf(command, "%s %s > %s\n", PDBCAT, pdb_name, PDBTMP);
    system(command); strcpy(file_name, PDBTMP);
  }else{
    sprintf(file_name, "%s", pdb_name);
  }
  FILE *file_in=fopen(file_name, "r");
  if(file_in==NULL){
    printf("\nWARNING, file %s not found\n", file_name); return(-1);
  }
  int N_model=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string, "MODEL", 5)==0)N_model++;
  }
  fclose(file_in);
  if(N_model==0)N_model=1;
  return(N_model);
}
/***************************************************************************/
/*                           Contact matrices                              */
/***************************************************************************/

struct contact *
Compute_cont_list(int *n_cont, int *Nc,
		  struct residue *res, int N_res, int ij_min)
{
  struct contact *cont_list=malloc(N_res*40*sizeof(struct contact)), *cont;

  if(init_map==0){
    cont_thr_c2=cont_thr_c*cont_thr_c;
    init_map++;
  }

  (*Nc)=0; cont=cont_list;
  short i_res, j_res;
  for(i_res=0; i_res<N_res; i_res++)n_cont[i_res]=0;
  for(i_res=0; i_res<N_res; i_res++){
    struct residue *res_i=res+i_res;
    for(j_res=i_res+ij_min; j_res< N_res; j_res++){
      if(Contact(res_i, res+j_res)){
	cont->res1=i_res; cont->res2=j_res; cont++;
	(*Nc)++; n_cont[i_res]++; 
      }
    }
  }
  struct contact *c_list=malloc((*Nc)*sizeof(struct contact)); int i;
  for(i=0; i<(*Nc); i++)c_list[i]=cont_list[i];
  free(cont_list);


  printf("Contact type MIN, thr= %.2f, %d native contacts\n",
	 cont_thr_c, *Nc);
  return(c_list);
}

struct contact *
Compute_cont_list_2(int *n_cont, int *Nc, float **min_dist,
		    struct residue *res, int N_res, int ij_min)
{
  struct contact *cont_list=malloc(N_res*40*sizeof(struct contact)), *cont;

  if(init_map==0){
    cont_thr_c2=cont_thr_c*cont_thr_c;
    init_map++;
  }

  (*Nc)=0; cont=cont_list;
  short i_res, j_res;
  for(i_res=0; i_res<N_res; i_res++)n_cont[i_res]=0;
  for(i_res=0; i_res<N_res; i_res++){
    struct residue *res_i=res+i_res; float *m=min_dist[i_res];
    for(j_res=i_res+ij_min; j_res< N_res; j_res++){
      m[j_res]=Min_dist(res_i, res+j_res);
      min_dist[j_res][i_res]=m[j_res];
      if(m[j_res]<=cont_thr_c){
	cont->res1=i_res; cont->res2=j_res; cont++;
	(*Nc)++; n_cont[i_res]++; 
      }
    }
  }
  struct contact *c_list=malloc((*Nc)*sizeof(struct contact)); int i;
  for(i=0; i<(*Nc); i++)c_list[i]=cont_list[i];
  free(cont_list);


  printf("Contact type MIN, thr= %.2f, %d native contacts\n",
	 cont_thr_c, *Nc);
  return(c_list);
}


int Contact(struct residue *res_i, struct residue *res_j){
  atom *atom1=res_i->atom_ptr, *atom2;
  float dx, dy, dz, *r1, *r2; int i, j;

  //printf("%c %c\n", res_i.amm, res_j.amm);
  for(i=0; i<res_i->n_atom; i++){
    r1=atom1->r;
    atom2=res_j->atom_ptr;
    for(j=0; j<res_j->n_atom; j++){
      r2=atom2->r;
      dx=(*(r1)  -*(r2));   if(fabs(dx)>cont_thr_c) goto new;
      dy=(*(r1+1)-*(r2+1)); if(fabs(dy)>cont_thr_c) goto new;
      dz=(*(r1+2)-*(r2+2)); if(fabs(dz)>cont_thr_c) goto new;
      if((dx*dx+dy*dy+dz*dz)<=cont_thr_c2) return(1);
    new:
      atom2++;
    }
    atom1++;
  }
  return(0);
}

float **Compute_min_dist(struct residue *res, int N_res)
{
  float **min_dist=malloc(N_res*sizeof(float *));
  for(short i_res=0; i_res<N_res; i_res++){
    min_dist[i_res]=malloc(N_res*sizeof(float));
  }
  for(short i_res=0; i_res<N_res; i_res++){
    float *m=min_dist[i_res];
    struct residue *res_i=res+i_res;
    for(short j_res=i_res+2; j_res< N_res; j_res++){
      m[j_res]=Min_dist(res_i, res+j_res);
      min_dist[j_res][i_res]=m[j_res];
    }
  }
  return(min_dist);
}

float Min_dist(struct residue *res_i, struct residue *res_j){
  atom *atom1=res_i->atom_ptr, *atom2;
  float dx, dy, dz, d2, *r1, *r2; int i, j;
  float min=-1;

  //printf("%c %c\n", res_i.amm, res_j.amm);
  for(i=0; i<res_i->n_atom; i++){
    r1=atom1->r;
    atom2=res_j->atom_ptr;
    for(j=0; j<res_j->n_atom; j++){
      r2=atom2->r;
      dx=(*(r1)  -*(r2));
      dy=(*(r1+1)-*(r2+1));
      dz=(*(r1+2)-*(r2+2));
      d2=dx*dx+dy*dy+dz*dz;
      if((min<0)||(d2<min))min=d2;
      atom2++;
    }
    atom1++;
  }
  if(min<0)return(100);
  return(sqrt(min));
}


int Count_proteins_seq(FILE *file_in){
  int n=0; char string[1000];
  fgets(string, sizeof(string), file_in);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')n++;
  }
  return(n);
}

int Get_monomeric(char *pdb_file){
  FILE *file_in=fopen(pdb_file, "r");
  if(file_in==NULL){
    printf("ERROR, PDB file %s does not exist\n", pdb_file);
    return(-1);
  }
  int monomeric=0; char string[200];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string, "REMARK 350", 10)==0){
      //REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC
      //REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DIMERIC
      char *s=string, *start, meric[80];
      while(*s!='\0'){
	if(*s==' '){
	  start=s+1;
	}else if(strncmp(s, "MERIC", 5)==0){
	  sscanf(start, "%s", meric);
	  if(strncmp(meric, "MONOMERIC", 9)==0)monomeric=1;
	  printf("%s oligomerization state: %s\n", pdb_file, meric);
	  break;
	}
	s++;
      }
    }else if(strncmp(string, "ATOM", 4)==0){
      break;
    }
  }
  fclose(file_in);
  return(monomeric);
}

int Get_NMR(char *pdb_file){
  FILE *file_in=fopen(pdb_file, "r");
  if(file_in==NULL){
    printf("ERROR, PDB file %s does not exist\n", pdb_file);
    return(-1);
  }
  int nmr=0; char string[200];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string, "EXPDTA", 6)==0){
      char *s=string;
      while(*s!='\0'){
	if(strncmp(s, "NMR", 3)==0){nmr=1; break;}
	s++;
      }
    }else if(strncmp(string, "REMARK", 6)==0){
      break;
    }
  }
  fclose(file_in);
  return(nmr);
}
