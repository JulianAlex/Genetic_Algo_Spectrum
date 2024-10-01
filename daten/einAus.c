#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/****************************************************************

Programm zum formatierten ausgeben

Aufruf: "einAus in.dat"

*****************************************************************/

#define ECHARGE 1.60217733*pow(10,-19)

main(int argc, char *argv[ ]){       

  FILE *in_file;
  char *cmd, *in_dat_1;
  int i, ch, n1=0;

  if(argc<2)
  {
    printf( "Missing argument\n" );
    return 1;
  }

  if((cmd = malloc(strlen(argv[1])+1)) == NULL)
  {
    printf( "Out of memory\n" );
    return 1;
  }


  in_dat_1 = argv[1];


  // === Zählen der Punkte ===========================

  in_file = fopen(in_dat_1, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat_1);
    exit(8);
  }
  while (1) {            //zaehle zeilenzahl = atomzahl
    ch = fgetc(in_file);
    if (ch == '\n')
      ++n1;
    if (ch == EOF)
      break;
  } 
  fclose(in_file);


  double x1[n1], x2[n1], x3[n1], x4[n1], x5[n1], x6[n1], x7[n1], x8[n1], x9[n1], x10[n1];
  

  // === Einlesen des Spektrums =========================

  in_file = fopen(in_dat_1, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat_1);
    exit(8);
  }
  for(i=0; i<n1; i++){
    //    fscanf(in_file, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
    //	   &x1[i], &x2[i], &x3[i], &x4[i], &x5[i], &x6[i], &x7[i], &x8[i], &x9[i], &x10[i]); 
    fscanf(in_file, " %lf %lf\n", &x1[i], &x2[i]); 
  }
  fclose(in_file);  


 // === Formatierte Ausgabe Spektrums =========================

  for(i=0; i<n1; i++)
    //   printf(" %10.3lf %14.10lf %14.10lf %14.10lf \n",  x1[i], x2[i], x3[i], x4[i]); 
    printf(" %10.3lf %14.10lf\n", x1[i], x2[i]); 

}
