#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>							
#include <fftw3.h>  

#include "eispack.h"     // Lineare Algebra, EW u EV berechnen 
#include "normal.h"      // Normalverteilte Zufallszahlen

#include "spek_funk.h"   // Funktions-Deklarationen


// Array
#define NTAU   101   //500     //101//  siehe "lineshapeArray( )" !!
#define DTAU   50.0  //10.0    //1.0//3.0//10.//6.// [fs] 
#define DTAU_0 0.001           //130.// 
#define NGAM   21    //81       //20   
#define DGAM   0.05  //0.0125   //0.025 
#define DGAM_0 0.10 //0.001//0.50 // gamma werte [0.5; 1.0] 

// Groebere Rasterung verfaelscht die Temperaturabhaengigkeit ein bisschen!!


// FFT-Parameter
#define NT 16384 //32768 //262144 //65536 //131072 //32768 // 2^14=16384    // anzahl element fuer FourTransf (2^n)!!
  
#define COUPLINGS 2 // 1 = calc Couplings in point dipole approx from coordinates
                    // 2 = read in couplings calculated elsewhere;

#define SPECDENS 1  // 1 = B777 spectral density
                  
#define SPECIES 1   // 1 = Tepidum-FMO    // geht in Spectraldichte ein, wenn SPECDEN > 1
                    // 2 = Aestuarii-FMO

#define BCL_8 0.35   // Anteil BCL 8, BCL7+BCL8 = 1

#define N 1000//000//000//      // N mal wuerfeln bei MonteCarlo, zum Rumprobieren reicht 10.000 => 5 Minuten
                                            // 100.000 => 30 min, e+6 => 5h , e+7 > 24h

#define TAU_0 1000.// 2750.     // [fs]   Pure Dephasing, ca. 5.3 cm-1 = 1000 fs, fuer Hole-Burning 2740.

#define LIFETIMEBROADENING 1    //1 Ein, 0 Aus

#define NNATM 4 // numb of N-Atom-Coord per pigment

//Grenzen für Spektrumsberechnung   
#define SPEKMIN 11800 //[cm-1]       // 11900-13000 => 2704
#define SPEKMAX 13000  
#define SPEKDIFF (SPEKMAX-SPEKMIN)

  
#define NEXC 200           // pigmExcDistrib()



// Float-Consts -----------------------------------------------------------------------------------------------

#define TEMP 4.0 //298.//77.      // [K]   Achtung: als float x.y schreiben!!  Tep 6K, Aest 4K

#define HRFAK 0.615 // Huang-Rhys-Skalierungsfaktor: so wählen, dass S1+S2 = Huang Rhys Faktor 

#define BETA 0.0    //Winkel zwischen qy und NB-ND (N-Atome im Pigment), nur PDIP(Punktdipol) !!!


//Vakuum Dipolstaerke fuer QY Uebergang [Einheit Debye] von 1 Pigment (µm aus denen dann µM berechnet wird)

#define DIPST sqrt(37.1) //3.6//4.0 // [DEBYE]   // Chl_a 4.0, Chl_b 3.6, BChl_a sqrt(37.1)
 


//-------------Site Energies----------------------------

// Mittelwerte Site-Energies
#define SITE_0 12500. //12445. // 12505. // 12450. // 12445. // [1/cm]
#define SITE_1 12425. //12520. // 12425. // 12450. // 12520.  
#define SITE_2 12195. //12205. // 12195. // 12450. // 12205.  
#define SITE_3 12375. //12335. // 12375. // 12450. // 12335.  
#define SITE_4 12600. //12490. // 12600. // 12450. // 12490.  
#define SITE_5 12515. //12640. // 12515. // 12450. // 12640.  
#define SITE_6 12465. //12450. // 12465. // 12450. // 12450.  
#define SITE_7 12700. //12800. // 12700. // 12450. // 12450. // (FMO-Tepidum new has 8 Chls)  

// 12421  12511  12196  12319  12509  12551  12472  12347 //Fit_2

// Breite der Gauss-Vert der site-Energies   // [1/cm]
#define FWHM_0 100. //0.01  // [1/cm]
#define FWHM_1 100. //0.01
#define FWHM_2 100. //0.01
#define FWHM_3 100. //0.01
#define FWHM_4 100. //0.01
#define FWHM_5 100. //0.01
#define FWHM_6 100. //0.01
#define FWHM_7 100. //0.01


// ------------ Burning ----------------------------------

#define E_BURN 12200. // [1/cm], Chl_b: 15047=664.6, 15239=656.2, 15293=653.9 // BurnFreq 
                      //         Chl_a: 14663=682.0, 14925=670.0, 15035=665.1   

#define WBZL 7.5// 2.5// 7.5 //10.//5.     // [1/cm] BurnWidth in ZL 

 
#define DIELEK 0.8  //Dielktr. Abschirm. Faktor, fuer PointDip-Calc und TransCharge (empty cavity faktor f)
                    //0.615 empty cav faktor fuer eps=2.56, 0.796 aus mflx mit eps=2 

#define NINDEX 1.414213562 // Brechungsindex, n = sqrt(epsilon), und epsilon=2

#define DIPFAK 1.0  //Skalieren der Dipolstärke unabh vom Dielektrikum

#define E_LAMB HRFAK*102 //[1/cm] Reorganisationsenergie E_lambda (4.17)

//--------------Fuer SpecDens J(w) aus Renger Marcus 2002 (2.16)

#define S1 HRFAK*0.8  // SpecDens Renger Marcus (2.16)
#define S2 HRFAK*0.5

#define W1 1.048294157e-4 // [fs-1]   // = 0.5565 1/cm * WZ2WFS => 1/fs Renger Marcus (2.16) 
#define W2 3.646240546e-4             // = 1.9357 1/cm * WZ2WFS => 1/fs Renger Marcus (2.16) 





#define W1CM 0.5565         // [cm-1] dasselbe wie oben
#define W2CM 1.9357         // [cm-1]

#define PI    3.141592654
#define TWOPI 6.283185307
#define PISQ  9.869604401  // PI^2

#define EXCLENG (1.0*(SPEKMAX-SPEKMIN))/NEXC


// FFT-Parameter
#define DELTA 1.0 //5.0 // t-steps FFT [fs]  // Freq-Aufloesung reziprok!  // Max. 5.0 !! 
#define WSTEP (TWOPI/(DELTA*NT))  // [fs-1] Schrittweite in Frequenz-Domaene



// Makros:

#define SQ(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define ABS(x) sqrt( ((x) * (x)) )



// Conversion Factors -----------------------------------------------------------------------------

#define EVWZ   8.06554                // (meV => 1/cm)
#define WZ2WFS 1.88369796e-4          // E(1/cm)*WZ2WFS = w(1/fs), WZ2WFS=2*Pi*c*100/e+15
#define ALPHFAC 4.340277778e-12  // Conversion [Debye^2*s/m]*4.34e-12 = e^2*Ang*s, 4.34e-12=10^-10/4.8^2



// Basic Constants

#define HBAR    0.6582119282           // eV*fs   (hbar/e)
#define HBDKB   7638.23301             // hquer/kB in fs*K
#define HBARC   1.973270524e-5         // hquer*c in meV*m 
#define CVAK    2.99792458e+8          // vacuum speed of light m/s

#define ECHARGE 1.60217733e-19
#define EPS_0   8.854187818e-12
#define ANG     1.0e-10
#define MILLI   0.001
#define DEBFAC  4.8         // Dipolstärke: e*Ang = 4.8 Debye

#define FWHMSIG 2.354820045 // 2*sqrt(2*ln(2)), fwhm=2*sqrt(2*ln(2))*sig (Gauss)  




