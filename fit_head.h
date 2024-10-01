#include "nr.h"
#include "nrutil.h"
#include "jutil.h"

#define NPOP 10//200//10//500//20     // Population size. ATTENTION: use 10, 20, ..., 50 etc. !!!
#define NFIT 3//100//3//100//100//  

// mit N=1000, NPOP=20, NFIT=100 dauerts 1.8 Stunden
// mit N=1000, NPOP=50, NFIT=100 dauerts 6.6 Stunden 

#define MUTMAX 30    
#define STARTMUTMAX 300 
#define REPMUTMAX 10
#define DIPMUTMAX 1

#define SP 1.2        // Selektionsdruck

// Mutationsrate soll ca 1% sein, bei 7 Eintraegen pro Chrom => 0.07
// (unterschied zu Reprod: wird nicht probabilistisch ausgewählt,
// sondern nach gewichtung) 
// Anteil der NICHT durch Nachkommen ersetzten Elterngeneration
// Summe ueber alle RATEs = 1  !!!

#define NINSRT_RATE 0.005      
#define REKOMB_RATE 0.01
#define REK_2_RATE  0.6
#define MUTAT_RATE  0.2      
#define REPROD_RATE (1-REKOMB_RATE-MUTAT_RATE-NINSRT_RATE)  


#define SITE_MIN 12100 
#define SITE_MAX 12800 

#define FWHM 100


#define FITMIN 12000 //11900 // Bereich, fuer den das Spektrum angepasst wird
#define FITMAX 12800 //13000

#define FITCUT 20 // Schneidet die ersten und letztens FITCUT elemente ab, 
                  // da am Rand Probleme mit der Ableitung

#define FAC_OD 3. // Wichtungsfaktor fuer OD-Spektrum
#define FAC_CD 1.
#define FAC_LD 1.
#define FAC_DX 2. // Wichtungsfaktor fuer OD'-Spektrum (Ableitung von OD)
