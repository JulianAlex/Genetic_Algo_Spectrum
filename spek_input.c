#include "spek_head.h"  
  

int countInput(char *in_dat){   
  
  // reads in input-file and counts number of input-lines

  FILE *in_file;

  int ch=0, n=0;

  in_file = fopen(in_dat, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(8);
  }
  while (1) {            
    ch = fgetc(in_file);
    if (ch == '\n')
      ++n;
    if (ch == EOF)
      break;
  }
  
  fclose(in_file);
  printf("input-lines: %d\n", n); printf("\n"); 

  return n;

}

 
// ==========================================================================

void readCoordinates(int *num, double *x, double *y, double *z, int n, char *in_dat){

  FILE *in_file;
  int i=0;  
 
 
  in_file = fopen(in_dat, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat);  
    exit(8); 
  }
  for(i=0; i<n; i++)
     fscanf(in_file, "%d %lf %lf %lf",
	    &num[i], &x[i], &y[i], &z[i]); 
  fclose(in_file);

}

// ==========================================================================

void readTRImerCoord(int *trim_num, double **trim_r, int trim_nz, char *in_dat){

  FILE *in_file;
  int i=0;
 
 
  in_file = fopen(in_dat, "r"); 
  if (in_file == NULL) { 
    (void)printf("Can not open %s\n", in_dat);
    exit(8);
  }
  for(i=0; i < trim_nz; i++){
    fscanf(in_file, "%d %lf %lf %lf\n",
	   &trim_num[i], &trim_r[i][0], &trim_r[i][1], &trim_r[i][2]); }
  fclose(in_file);

}

// ==========================================================================

void readDipoles(int nz, double **qy, char *in_dat_1, char *in_dat_2){

  FILE *in_file; 

  in_file = fopen(in_dat_1, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat_1);  
    exit(8); 
  }
  fscanf(in_file, "%lf %lf %lf", &qy[0][0], &qy[0][1], &qy[0][2]); 
  fclose(in_file);


   
  in_file = fopen(in_dat_2, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat_2);  
    exit(8); 
  }
  fscanf(in_file, "%lf %lf %lf", &qy[1][0], &qy[1][1], &qy[1][2]);   
  fclose(in_file);

  printf("Dipole moments from files:  %s,  %s\n\n", in_dat_1, in_dat_2);

}

// ==========================================================================

void readDipol2(int nz, double **qy, char *in_dat_1){

  FILE *in_file; 
  int i; 

  in_file = fopen(in_dat_1, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat_1);  
    exit(8); 
  }
  for(i=0; i<nz; i++)
  {
    fscanf(in_file, "%lf %lf %lf", &qy[i][0], &qy[i][1], &qy[i][2]); 
  }
  fclose(in_file);

  printf("Dipole moments from file:  %s\n\n", in_dat_1);

}
//============================================================================

void readCouplings(int nz, double **vab){

  int i=0, j=1;

  double v_12 = 0; //COUP_12*DIPFAK;

  printf("Couplings from files: spek.h \n\n");

  vab[0][1] = v_12;

  for(i=0; i<nz; i++)
    vab[i][i]=0.0;

  for(i=1; i<nz; i++)
    for(j=0; j<i; j++)
      vab[i][j] = vab[j][i];

  for(i=0; i<nz; i++)
    for(j=0; j<nz; j++)
      vab[i][j] *= EVWZ;

  printf("Couplings in cm-1:\n");
  for(i=0; i<nz; i++){
    for(j=0; j<nz; j++){
      printf(" %11.4f ", vab[i][j]);
    }printf("\n");
  }printf("\n");


  printf("Couplings in meV:\n");
  for(i=0; i<nz; i++){
    for(j=0; j<nz; j++){
      printf(" %11.4f ", vab[i][j]/EVWZ);
    }printf("\n");
  }printf("\n");


}

//================================================================================================

// function "eingabeExpSpek" reads in experimental spektra from file 

 
void eingabeExpSpek_OD(char *in_dat, double *exp_ab, int nx){

  FILE *in_file;

  int i;
  double w;
 
  in_file = fopen(in_dat, "r");
  if (in_file == NULL) {
     (void)printf("Can not open %s\n", in_dat);
      exit(8);
  }

  for(i=0; i<nx; i++)
  {
    fscanf(in_file, "%lf %lf\n", &w, &exp_ab[i]);  
  }
  fclose(in_file);

  
  printf("EingabeFile:  %s\n\n", in_dat);
}


//===============================================================================================

// function "eingabeExpSpek" reads in experimental spektra from file 

 
void eingabeExpSpek(char *in_dat, double *exp_ab, double *exp_cd, double *exp_ld, double *exp_dx, int nx){

  FILE *in_file;

  int i;
  double w;
 
  in_file = fopen(in_dat, "r");
  if (in_file == NULL) {
     (void)printf("Can not open %s\n", in_dat);
      exit(8);
  }

  for(i=0; i<nx; i++)
  {
    fscanf(in_file, "%lf %lf %lf %lf %lf\n", &w, &exp_ab[i], &exp_cd[i], &exp_ld[i], &exp_dx[i]);  
  }
  fclose(in_file);

  
  printf("EingabeFile:  %s\n\n", in_dat);
}

