// Julian Adolphs
// File "fitnVert.c"

// void diffExpSimVSB(...)
// void verteilungsfkt(...)

#include "spek_head.h"
#include "fit_head.h"
 
//===========================================================================

void diffExpSim(int i, int anz, double wstep, double *fitness, double *nr_fitness, 
		double *sumspek_ab, double *sumspek_cd, double *sumspek_ld, double *sumspek_dx, 
		double *exp_ab, double *exp_cd, double *exp_ld, double *exp_dx){
    
  int k, kmin, kmax; 
  double diff = 0.0, diff_1 = 0.0, diff_2 = 0.0, diff_3 = 0.0, diff_4 = 0.0;

  kmin = (int) ((FITMIN-SPEKMIN)/wstep + 0.5 + FITCUT); // an den Raendern gibts Probleme....
  kmax = (int) ((FITMAX-SPEKMIN)/wstep + 0.5 - FITCUT);   

  //  printf("(fit_fitnessFkt): k, kmax, wstep %d  %d  %lf\n", k, kmax, wstep);

  for(k = kmin; k < kmax; k++) 
  {
    diff_1 += sqrt( SQ( exp_ab[k]  - sumspek_ab[k] ) );
    diff_2 += sqrt( SQ( exp_cd[k]  - sumspek_cd[k] ) );
    diff_3 += sqrt( SQ( exp_ld[k]  - sumspek_ld[k] ) );
    diff_4 += sqrt( SQ( exp_dx[k]  - sumspek_dx[k] ) );
  }     


  diff = FAC_OD*diff_1 + FAC_CD*diff_2 + FAC_LD*diff_3 + FAC_DX*diff_4;  //FAC_... Wichtungsfaktoren

  fitness[i] = 1.0/diff;  

  nr_fitness[i+1] = diff;

  // printf("diffExpSim:  %d  %g %g %g %g\n", i, FAC_OD*diff_1, FAC_CD*diff_2, FAC_LD*diff_3, FAC_DX*diff_4);


} 

//===========================================================================

void verteilungsfkt(double *fitness, double *nr_fitness, double *rank_fitn, double *vertlgsfkt, 
		    double *prob, int *irank, int *nr_irank, double **parent){

// PROCEDURE "verteilungsfkt" berechnet eine Verteilungsfunktion nach der
// Methode "Lineares Ranking" nach W. Kinnebrock: 
// Optimierung mit genetischen und selektiven Algorithmen (S.73)
// und macht die Ausgabe von Rang, Fitnesswert und Nachkommen in eine
// Datei. (evtl noch durch eine andere Procedure erledigen lassen!)

  FILE *out_file;
  char *out_dat="fitAes78_1.out";
  out_file = fopen(out_dat, "a"); 

  int j, m;  
  double sum = 0.0, max = 1.99, min = 0.01, n = NPOP, probsum = 0.0;

  for(m=0; m<NPOP; m++)
    irank[m] = nr_irank[m+1];

  for(m=0; m<NPOP; m++)
    rank_fitn[m] = 1/nr_fitness[1+m];
 
  for(m=0; m<NPOP; m++)
    sum += fitness[m];

  fprintf(out_file, "MaxFitn = %g, MinFitn = %g, GesFitn = %g\n", 
	  rank_fitn[0], rank_fitn[NPOP-1], sum); fprintf(out_file, "\n");
 
  for(m=0; m<NPOP; m++)
      fprintf(out_file, "%3d %3d  %12.6g   %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n",
	      m, nr_irank[m+1], fitness[m], 
	      parent[m][0], parent[m][1], parent[m][2], parent[m][3], parent[m][4], parent[m][5], parent[m][6], parent[m][7]); 

  prob[0] = 0.0;   
  for(j=1; j<=NPOP; j++){  
    prob[j] = ( max - (max-min)*(j-1)/(n-1) )/n;  
  }
  
  for(j=1; j<=NPOP; j++)
    probsum += prob[j]; 

  sum = 0.0; 
  for(j=0; j<=NPOP; j++){  // Verteilungsfunktion geht durch den ursprung! 
    sum += prob[j];
    vertlgsfkt[j] = sum/probsum;   
  }

  fprintf(out_file, "\n");
  fclose(out_file); 
}

