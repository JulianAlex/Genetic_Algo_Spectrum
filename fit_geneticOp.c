// Julian Adolphs
// File "genOperation.c"
// Enthaelt die FUNKTIONEN und PROCEDUREN

// double mutElement(...)
// int selectElement(...)
// void selectNinsrt(...)
// void notinsert(...)
// void twoSelect(...)
// void rekomb2(...)
// void selection(...)
// void mutation(...)
// void rekombination(...)
// void reproduktion(...)
// void repromut(...)
// void geneticOperate(...)


#include "spek_head.h"
#include "fit_head.h"
 



//==============================================================================

double mutElement(int mut_elem, int nz, int i, double **offspring, 
		  double mutmax){

// Function "mutElement" berechnet um wieviel ein einzelnes Element mutiert
// werden soll. 

  double help_double, off_spring; 

  do
  {
    help_double = offspring[i][mut_elem] + mutmax*( 1 - 2*((double) rand()/RAND_MAX) );
  }
  while( (help_double < SITE_MIN) || (help_double > SITE_MAX) );   
 
  off_spring = help_double;

  return off_spring;
}


//---------------------------------------------------------------------------

int selectElement(int nz){

// Function "selectElement" waehlt probabilistisch ein Element zum
// Mutieren aus.

  int j, mutelem = 0; 
  double xj, pmut;

  pmut = (double) rand()/RAND_MAX;
 
  for(j=0; j<nz; j++)
  {
    xj = (double) j;  
    if( ( pmut >  xj/(nz) ) && ( pmut < (xj+1)/(nz) ) ) 
      mutelem = j;
  }
  return mutelem;  
}

//===========================================================================

 
void selectNinsrt(int *irank, int i, int k, int *num_sel ){

  // Der selektierte Index bezieht sich auf das ranking!!
  // muss jetzt wieder dem soundsovieltem Chromosom zugeordenet werden: 

  int j;  

  for(j=0; j<NPOP; j++)
    if(irank[j] == i+1)
    {
      *num_sel = j;
      break;
    }

}

//---------------------------------------------------------------------------

void notinsert(double *fitness,  int *irank, int i, int k, 
		  int nz, int *num_sel,
		  double **parent, double **offspring){

// Ersetzt einen Nachkommen durch einen Eltern (Klon :-)

  int j;  

  for(j=0; j<nz; j++)
    offspring[i][j] = parent[ *num_sel][j]; 
}
  
 
//----------------------------------------------------------------------------


void twoSelect(double *vertlgsfkt, int *irank, int i, int k, int *indx_sel_1,
	       int *num_sel_1, int *indx_sel_2, int *num_sel_2){

  double ran1, ran2;

  int j;  
 
  ran1 = (double) rand()/RAND_MAX;   // zufallszahl zwischen 0 und 1
  ran2 = (double) rand()/RAND_MAX; 

  for(j=0; j<NPOP; j++)
    if ( (vertlgsfkt[j]<= ran1) && (vertlgsfkt[j+1]>ran1) )
    {
      *indx_sel_1 = j+1;
    }

  for(j=0; j<NPOP; j++)
    if ( (vertlgsfkt[j]<= ran2) && (vertlgsfkt[j+1]>ran2) )
    {
      *indx_sel_2 = j+1;
    }

  // Der selektierte Index bezieht sich auf das ranking!!
  // muss jetzt wieder dem soundsovieltem Chromosom zugeordenet werden:

  for(j=0; j<NPOP; j++)
    if(irank[j] == *indx_sel_1)
    {
      *num_sel_1 = j;
    }

   for(j=0; j<NPOP; j++)
    if(irank[j] == *indx_sel_2)
    {
      *num_sel_2 = j;
    }

}

//---------------------------------------------------------------------------


void rekomb2(double *fitness, double *vertlgsfkt, int *irank, int i, int k, 
	     int nz, int *indx_sel_1, int *num_sel_1, int *indx_sel_2, 
	     int *num_sel_2, double **parent, double **offspring){

  int j, cut_elem=1, mutelem_1, mutelem_2;
 
  double pcut, xj,  repmutmax=REPMUTMAX;

  pcut = (double) rand()/RAND_MAX;

  // Es soll mind. nach dem ersten und hoechstens nach dem vorletzten element 
  // "cut and changed" werden:

  for(j=1; j<nz; j++)
  { 
    xj = (double) j;   // Type-cast
    if( ( pcut >  (xj-1)/(nz-1) ) && ( pcut < xj/(nz-1) ) )             
      cut_elem = j;
  } 
   
  for(j=0; j<cut_elem; j++)
  {
    offspring[i][j]   = parent[*num_sel_1][j];
    offspring[i+1][j] = parent[*num_sel_2][j];
  }
  for(j=cut_elem; j<nz; j++)
  {
    offspring[i][j]   = parent[*num_sel_2][j];
    offspring[i+1][j] = parent[*num_sel_1][j];
  }
   
 // Jetzt noch jeweils ein element mutieren: 
  
  mutelem_1 = selectElement(nz);
  mutelem_2 = selectElement(nz);
 
  offspring[i][mutelem_1] = mutElement(mutelem_1, nz, i,offspring,repmutmax);
  offspring[i+1][mutelem_2]=mutElement(mutelem_2,nz,i+1,offspring,repmutmax);

  /******
  printf(" Rekom_2_Mut:  %2d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n\n", 
	  mutelem_1, offspring[i][0], offspring[i][1], offspring[i][2], 
	  offspring[i][3], offspring[i][4], offspring[i][5], offspring[i][6]);

  printf(" Rekom_2_Mut:  %2d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n\n", 
	  mutelem_2, offspring[i+1][0], offspring[i+1][1], offspring[i+1][2], 
	 offspring[i+1][3], offspring[i+1][4], offspring[i+1][5], 
	 offspring[i+1][6]);
  *******/ 
  
}

//---------------------------------------------------------------------------


void selection(double *vertlgsfkt, int *irank, int i, int k, int *indx_sel, 
	       int *num_sel ){

  double ran;

  int j;  

  ran = (double) rand()/RAND_MAX;   // zufallszahl zwischen 0 und 1

  for(j=0; j<NPOP; j++)
    if ( (vertlgsfkt[j]<= ran) && (vertlgsfkt[j+1]>ran) )
    {
      *indx_sel = j+1;
      break;
    }

  //  printf("Selekt Chrom (nach ranking): %3d  %6.3f \n", *indx_sel, ran); 

  // Der selektierte Index bezieht sich auf das ranking!!
  // muss jetzt wieder dem soundsovieltem Chromosom zugeordenet werden:

  for(j=0; j<NPOP; j++)
    if(irank[j] == *indx_sel)
    {
      *num_sel = j;
      break;
    }

}

//---------------------------------------------------------------------------


void mutation(double *fitness, double *vertlgsfkt, int *irank, int i, int k, 
		  int nz, int *indx_sel, int *num_sel,
		  double **parent, double **offspring){

  double panz, mutmax=MUTMAX;
        
  int j,  mut_elem1 = 0, mut_elem2 = 0;

  panz = (double) rand()/RAND_MAX;
   
  for(j=0; j<nz; j++)
    offspring[i][j] = parent[*num_sel][j];   //selekt. Chrom. 
 
  if (panz < 0.5)                    //Mutation: waehle probab. EIN Elemente des Chromosoms...
  {                                  
    mut_elem1 = selectElement(nz);

    //...und veraendere es um einen zufaelligen Wert 

    offspring[i][mut_elem1] = mutElement(mut_elem1, nz, i, offspring, mutmax);
  }

  if (panz >= 0.5)                   //Mutation: waehle per zufall ZWEI Elemente des Chromosoms
  { 
    mut_elem1 = selectElement(nz);
    mut_elem2 = selectElement(nz);

    //...und veraendere es um einen zufaelligen Wert

    offspring[i][mut_elem1] = mutElement(mut_elem1, nz, i, offspring, mutmax);
    offspring[i][mut_elem2] = mutElement(mut_elem2, nz, i, offspring, mutmax);
  }

} 

//---------------------------------------------------------------------------


void rekombination(double *fitness, double *vertlgsfkt, int *irank, int i, 
		   int k, int nz, int *indx_sel_1, int *num_sel_1, 
		   int *indx_sel_2, int *num_sel_2,
		   double **parent, double **offspring){

  int j, mut_elem=0;
  double repmutmax=REPMUTMAX;

  //Intermediaere Rekombination (Mittelwertbildung) obiger beider chromosomen:

  for(j=0; j<nz; j++)
    offspring[i][j] = 0.5 * (parent[*num_sel_1][j] + parent[*num_sel_2][j]);

  // Jetzt noch ein element mutieren:

  mut_elem = selectElement(nz);

  offspring[i][mut_elem] = mutElement(mut_elem, nz, i, offspring, repmutmax);
 

  /********
  printf("ReKombMut: %d  %6.1f %6.1f  %6.1f  %6.1f  %6.1f  %6.1f  %6.3f\n\n", 
	 mut_elem, offspring[i][0],offspring[i][1],offspring[i][2],
	 offspring[i][3],offspring[i][4],offspring[i][5],offspring[i][6]);
  ********/
}


//---------------------------------------------------------------------------


void reproduktion(double *fitness, double *vertlgsfkt, int *irank, int i, 
		  int k, int nz, int *indx_sel, int *num_sel,
		  double **parent, double **offspring){

  int j;  

  for(j=0; j<nz; j++)
    offspring[i][j] = parent[ *num_sel][j]; 

  /********
  printf("Reprodukt: %d %6.3f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f\n\n", 
	 *num_sel, fitness[ *num_sel], 
	 offspring[i][0],offspring[i][1],offspring[i][2],offspring[i][3],
	 offspring[i][4],offspring[i][5],offspring[i][6]);
  ********/
}

//----------------------------------------------------------------------------


void repromut(double *fitness, double *vertlgsfkt, int *irank, int i, int k, 
		  int nz, int *indx_sel, int *num_sel,
		  double **parent, double **offspring){
 
  double panz, mutjanein, repmutmax=REPMUTMAX;
        
  int j,  mut_elem1 = 0, mut_elem2 = 0;

  panz = (double) rand()/RAND_MAX;
  mutjanein = (double) rand()/RAND_MAX;

  for(j=0; j<nz; j++)
    offspring[i][j] = parent[*num_sel][j];   //selekt. Chrom. 
 
  
 
  if (mutjanein < 0.2)
  {
    if (panz < 0.5)    //Mutation: waehle per zufall EIN Elemente des Chromosoms
    {  
      mut_elem1 = selectElement(nz);

      // ...und veraendere es um einen zufaelligen Wert zwischen 
      // -REPMUTMAX und +REPMUTMAX !!!
      // unter der bedingung dass es aus [SITE_MIN, SITE_MAX] ist

      offspring[i][mut_elem1] = mutElement(mut_elem1, nz, i, offspring, 
					   repmutmax);
      /********
      printf("ReproMut_1: %d %d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f\n\n",
	     *num_sel, mut_elem1, 
	     offspring[i][0],offspring[i][1],offspring[i][2],offspring[i][3],
	     offspring[i][4],offspring[i][5],offspring[i][6]);
      ********/
    }
    

    if (panz >= 0.5){ //Mutation: waehle per zuf ZWEI Elem des Chromosoms...

      mut_elem1 = selectElement(nz);
      mut_elem2 = selectElement(nz);

      //...und veraendere es um einen zufaelligen Wert zwischen 
      //-REPMUTMAX und +REPMUTMAX !!!
      // unter der bedingung dass es aus [SITE_MIN, SITE_MAX] ist


      offspring[i][mut_elem1] = mutElement(mut_elem1, nz, i, offspring, 
					   repmutmax);
      offspring[i][mut_elem2] = mutElement(mut_elem2, nz, i, offspring, 
					   repmutmax);

    
      /********
     printf("RepMut_2: %d %d %d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f\n\n",
	     *num_sel, mut_elem1, mut_elem2, 
	     offspring[i][0],offspring[i][1],offspring[i][2],offspring[i][3],
	     offspring[i][4],offspring[i][5],offspring[i][6]);
      ********/
    }
  }  
}

//=================================================================================== 


void geneticOperate(int n1, int n2, int n3, int n4, int nz, int i, int k, 
		    int *irank, int *num_sel, int *num_sel_1, int *num_sel_2, 
		    int *indx_sel, int *indx_sel_1, int *indx_sel_2,
		    double *fitness, double *vertlgsfkt, double **parent, double **offspring){


  for( i = 0;  i < n1;  i++ )
  {
    selectNinsrt(irank, i, k, num_sel); 
    notinsert(fitness, irank, i, k, nz, num_sel, parent, offspring);  //klonen 
  } 
  for(i = n1;  i < n2;  i++)
  {   
    twoSelect(vertlgsfkt,irank,i,k,indx_sel_1,num_sel_1,indx_sel_2,num_sel_2);
    rekombination(fitness, vertlgsfkt, irank, i, k, nz, indx_sel_1, num_sel_1, indx_sel_2, num_sel_2, parent, offspring); 
  }         
  for(i = n2;  i < n3;  i=i+2)
  {   
    twoSelect(vertlgsfkt,irank,i,k,indx_sel_1,num_sel_1,indx_sel_2,num_sel_2);
    rekomb2(fitness, vertlgsfkt, irank, i, k, nz, indx_sel_1, num_sel_1, indx_sel_2, num_sel_2, parent, offspring); 
  } 
  for(i = n3;  i < n4;  i++)
  {
    selection(vertlgsfkt, irank, i, k, indx_sel, num_sel);
    mutation(fitness, vertlgsfkt, irank, i, k, nz, indx_sel, num_sel, parent, offspring);
  }    
  for( i = n4;  i < NPOP;  i++ )
  {
    selection(vertlgsfkt, irank, i, k, indx_sel, num_sel);
    reproduktion(fitness, vertlgsfkt, irank, i, k, nz, indx_sel, num_sel, parent, offspring);
    repromut(fitness, vertlgsfkt, irank, i, k, nz, indx_sel, num_sel, parent, offspring);
  }

  printf(" \n"); 
   
  int j, m;
  for(j=0; j<NPOP; j++)  
    for(m=0; m<nz; m++)
      parent[j][m] = offspring[j][m];  
  
}
