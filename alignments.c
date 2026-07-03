#include "alignments.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "NeedlemanWunsch.h"

static int Check_MSA(int *L, char *file_ali);

char **Read_MSA(int *N, int *L, char ***name_seq, int **selected,
		char *file_ali, float thr)
{
  if(file_ali==NULL)return(NULL);
  *N=Check_MSA(L, file_ali);
  char **MSA=malloc(*N*sizeof(char *));
  *name_seq= malloc(*N*sizeof(char *));
  int NCHAR=100, n;
  for(n=0; n<*N; n++){
    MSA[n]=malloc(*L*sizeof(char));
    (*name_seq)[n]=malloc(NCHAR*sizeof(char));
  }

  FILE *file_in=fopen(file_ali, "r");
  char string[10000], *seq=NULL; n=-1;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      n++; seq=MSA[n];
      sscanf(string+1, "%s", (*name_seq)[n]);
    }else{
      char *s=string;
      while((*s!='\n')&&(*s!='\0')){*seq=*s; seq++; s++;}
    }
  }
  fclose(file_in);

  *selected=malloc(*N*sizeof(int)); int m=0;
  for(n=0; n<*N; n++){
    int l=0, k; for(k=0; k<*L; k++)if(MSA[n][k]!='-')l++;
    if(l > (*L)*thr){(*selected)[n]=1; m++;}
    else{(*selected)[n]=0;}
  }

  printf("%d sequences with %d letters read in file %s, %d selected\n",
	 *N,*L,file_ali, m);
  return(MSA);
}

int Check_MSA(int *L, char *file_ali)
{
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("ERROR, MSA file %s does not exist\n", file_ali); exit(8);
  }

  printf("Reading MSA file %s\n", file_ali);
  char string[10000]; int n=-1, l=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      if(n==0){
	*L=l;
      }else if(l!=*L){
	printf("ERROR, different sequence lengths: 0 %d %d %d\n", *L, n, l);
	exit(8);
      }
      n++; l=0;
    }else{
      char *s=string;
      while((*s!='\n')&&(*s!='\0')){l++; s++;}
    }
  }
  if(n==0){
    *L=l;
  }else if(l!=*L){
    printf("ERROR, different sequence lengths: %d %d\n", *L, l);
    exit(8);
  }
  n++;
  fclose(file_in);
  if(n==0){
    printf("ERROR, no sequence found in MSA file %s\n", file_ali);
    exit(8);
  }
  return(n);
}

