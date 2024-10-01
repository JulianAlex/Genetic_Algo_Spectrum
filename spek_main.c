/****************************************************************************

Julian Adolphs

(Contributions of Manuel Berrer to hole burning)

Die inhomogenen Spektren in dynamischer Theorie werden mit einer MonteCarlo-Methode berechnet.

Der Input besteht aus den X-Ray-Koordinaten 
und den Site-Energies und den Dipolstaerken (spek.h)

Output: out_linSpectra.dat

Files die zum Programm "spektrum" gehoeren:  spek_*.c spek_*.h

Fourier-Transformationen mit FFTW-Routinen 

Diagonalisierung der Hamilton-Matrix mit "Eispack": eispack.c eispack.h

Zufallszahlen: normal.h normal.c

****************************************************************************/
  

#include "spek_head.h"   // Makros & Definitionen



int main(int argc, char *argv[ ]) 
{

  clock_t prgstart, prgende;  // Laufzeitmessung
  prgstart = clock();         // CPU-Zeit zu Beginn des Programmes

  int i=0, k=0, n=0, nz=0, hn=0, nrow=0, ncol=0, pigm_no=0; 

  int h_var = 17;        // seed for random number generator
  int *seed  = &h_var;

  int couplings = COUPLINGS;   // Kopplungen berechnen oder einlesen?!
  int dipoles   = DIPOLES;     // qy-Dipol-Richtung berechnen oder einlesen?!
  int wscp      = WSCP;        // Chl_ or Chl_b-WSCP?!  

  char *in_dat = "infile_koord.dat";
  n = countInput(in_dat);

  nz = n/NNATM;                  // anzahl der chlorophylle, WSCP:2
  hn = SQ(nz);
  printf("# Anzahl N-Atome:  %d ,  Anzahl Chl:  %d,  Temperatur:  %4.1lf K\n", n, nz, TEMP);
  
  double wstep = TWOPI/(NT*DELTA*WZ2WFS);   // [cm-1]     
  int    anz   = (int) (SPEKDIFF/wstep + 1);    
  printf("# wstep:  %lf cm-1, anz:  %d\n", wstep, anz);


 // ------------- Declare & Allocate vectors  ---------------------------------

  int   *num = calloc( n, sizeof(int));

  double *x = calloc( n, sizeof(double));
  double *y = calloc( n, sizeof(double));
  double *z = calloc( n, sizeof(double));

  double *dipst   = calloc( nz, sizeof(double));
  double *site    = calloc( nz, sizeof(double));   // Site Energie Mittelwerte 
  double *sigma   = calloc( nz, sizeof(double));
  double *alpha   = calloc( nz, sizeof(double));   // alpha=|mu|^2 dipst^2 [in Debye^2] 
  double *rotst   = calloc( nz, sizeof(double));   // rotational strength for CD 
  double *tau     = calloc( nz, sizeof(double));
  double *w_m0_sh = calloc( nz, sizeof(double));
  double *w_hb_sh = calloc( nz, sizeof(double));

  double *jw    = calloc( NT/2, sizeof(double));   // J(w)   jw[] = jw[0] ... jw[NT/2-1]

  double *cw_re = calloc( NT, sizeof(double));
  double *cw_im = calloc( NT, sizeof(double));

  double *spek_ab       = calloc( anz+1, sizeof(double)); // Homogenes Absorptions-Spektrum 
  double *spek_sb       = calloc( anz+1, sizeof(double)); // Homogenes Absorptions-Spektrum nur der SeitenBande
  double *spek_zl       = calloc( anz+1, sizeof(double)); // Homogenes Absorptions-Spektrum nur der ZerozeroLine
  double *spek_cd       = calloc( anz+1, sizeof(double)); // Homogenes CD-Spektrum 
  double *spek_hb       = calloc( anz+1, sizeof(double)); // HoleBurning Spektrum nach Brennvorgang an 1 Pigment
  double *spek_hb_res   = calloc( anz+1, sizeof(double)); // Resonant-HoleBurning Spektrum nach Brennvorgang an 1 Pigment
  double *sumspek_ab    = calloc( anz+2, sizeof(double)); // Inhomogenes Absorptionsspektrum (=aufsummiert)
  double *sumspek_sb    = calloc( anz+2, sizeof(double)); // Inhomogenes AbsSpek der SB (=aufsummiert)
  double *sumspek_zl    = calloc( anz+2, sizeof(double)); // Inhomogenes AbsSpek der ZL (=aufsummiert)
  double *sumspek_cd    = calloc( anz+2, sizeof(double)); // Inhomogenes CD-spektrum (=aufsummiert)
  double *sumspek_hb    = calloc( anz+2, sizeof(double)); // Inhomog HoleBurning Spektrum 
  double *sumspek_hbab  = calloc( anz+2, sizeof(double)); // Skaliertes HoleBurning-Absorptionsspektrum
  double *sumspek_hbdif = calloc( anz+2, sizeof(double)); // Skaliertes HoleBurning-Absorptionsspektrum

  // Achtung: fftw_malloc initialisiert den Speicher nicht! 
  fftw_complex *gt; gt = (fftw_complex*) fftw_malloc( NT * sizeof(fftw_complex) );   // G(t)   
  fftw_complex *ct; ct = (fftw_complex*) fftw_malloc( NT * sizeof(fftw_complex) );   // C(t)   
  fftw_complex *cw; cw = (fftw_complex*) fftw_malloc( NT * sizeof(fftw_complex) );   //~C(w) 

  for(i=0; i<NT; i++)
  {
    gt[i] = 0.0;
    ct[i] = 0.0;
    cw[i] = 0.0;
  }

  fftw_complex *in, *out;
  fftw_plan plan_bwd, plan_fwd;

  in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NT); //allocating memory
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NT);

  for(k=0; k<NT; k++)
  {
    in[k]  = 0.0;
    out[k] = 0.0;
  }

  //  plan = fftw_plan_dft_1d(NT, in, out, +1, FFTW_MEASURE);  // -1/+1 direction of FT, FFTW_MEASURE accuracy
  plan_bwd = fftw_plan_dft_1d(NT, in, out, +1, FFTW_MEASURE);  // backward-transform
  plan_fwd = fftw_plan_dft_1d(NT, in, out, -1, FFTW_MEASURE);  // forward-transform


  //---- Matrizen -----------------------------------------------------------------------------------
  // Im FFTW-Manual wird von dieser Art Arrays zu konstruieren abgeraten, Schmaranz empfiehlt es so...

  double **mg = calloc( nz, sizeof(double*) );
  double **na = calloc( nz, sizeof(double*) );
  double **nb = calloc( nz, sizeof(double*) );
  double **nc = calloc( nz, sizeof(double*) );
  double **nd = calloc( nz, sizeof(double*) );
  double **qx  = calloc( nz, sizeof(double*) );
  double **qy   = calloc( nz, sizeof(double*) );
  double **vab   = calloc( nz, sizeof(double*) );
  double **vab_hb = calloc( nz, sizeof(double*) );
  double **rij   = calloc( nz, sizeof(double*) );
  double **gam  = calloc( nz, sizeof(double*) );
  double **j_w = calloc( nz, sizeof(double*) );     // Spectral Density J(w)
  double **n_w = calloc( nz, sizeof(double*) );     // No of vib.quanta n(w)

  nrow = nz;
  while(nrow--)
  {
    mg[nrow] = calloc(3, sizeof(double));
    na[nrow] = calloc(3, sizeof(double));
    nb[nrow] = calloc(3, sizeof(double));
    nc[nrow] = calloc(3, sizeof(double));
    nd[nrow] = calloc(3, sizeof(double));
    qx[nrow] = calloc(3, sizeof(double));
    qy[nrow] = calloc(3, sizeof(double));
    vab[nrow] = calloc(nz, sizeof(double));
    vab_hb[nrow] = calloc(nz, sizeof(double));
    rij[nrow] = calloc(nz, sizeof(double));
    gam[nrow] = calloc(nz, sizeof(double));
    j_w[nrow] = calloc(nz, sizeof(double));
    n_w[nrow] = calloc(nz, sizeof(double));
  }

  double **npig = calloc( nz, sizeof(double*) );
  double **nexc = calloc( nz, sizeof(double*) );
  nrow = nz;
  while(nrow--)
  {
    npig[nrow]  = calloc( NEXC, sizeof(double));
    nexc[nrow]  = calloc( NEXC, sizeof(double));
  }

  double **delta = calloc( hn, sizeof(double*));
  nrow = hn;
  while(nrow--)
    delta[nrow] = calloc(3, sizeof(double));
  
  double **arrayspek_sb = calloc( nz, sizeof(double*) );
  double **arrayspek_zl = calloc( nz, sizeof(double*) );

  nrow = nz;
  while(nrow--)
  {
    arrayspek_sb[nrow] = calloc( 4*anz, sizeof(double) ); 
    arrayspek_zl[nrow] = calloc( 4*anz, sizeof(double) );
  }

 
  // float statt double halbiert den benötigten Arbeitsspeicher
  float ***sideband = calloc( NGAM, sizeof(float**) );          //[gamma][tau][w]
  float ***zeroline = calloc( NGAM, sizeof(float**) );          //[gamma][tau][w]
  nrow = NGAM;
  while(nrow--){
    sideband[nrow] = calloc( NTAU, sizeof(float*) );
    zeroline[nrow] = calloc( NTAU, sizeof(float*) );
    ncol = NTAU;
    while(ncol--){
      sideband[nrow][ncol] = calloc( NT, sizeof(float));
      zeroline[nrow][ncol] = calloc( NT, sizeof(float)); 
    }
  }


  
  double *hamilt    = calloc( nz*nz, sizeof(double));    // Hamilton-Matrix
  double *ev_h      = calloc( nz*nz, sizeof(double));    // Eigenvektoren (1D) output of the eispack-routines
  double *eigval    = calloc( nz, sizeof(double));       // Eigenwerte    
  double *eigval_hb = calloc( nz, sizeof(double));       // Eigenwerte    

  double **eigvec    = calloc(nz, sizeof(double*));    // Eigenvektoren (2D)
  double **eigvec_hb = calloc(nz, sizeof(double*)); 
  nrow = nz;
  while(nrow--)
  {
    eigvec[nrow]    = calloc(nz, sizeof(double));
    eigvec_hb[nrow] = calloc(nz, sizeof(double));
  }




 //======= Begin Hauptprogramm "spektrum" ====================================

  in_dat="infile_koord.dat";
  readCoordinates(num, x, y, z, n, in_dat);
  parameter(site, sigma, dipst, nz);        
  koord(nz, x, y, z, mg, na, nb, nc, nd);
  pigmCenters(nz, mg, nb, nd, delta, rij);   

  if( (couplings == 1) && (dipoles == 1) )
  {
    dipolmom(nz, na, nb, nc, nd, qx, qy);
    interaction(nz, qy, delta, vab, rij, dipst);
    printf("PointDipoleApproximation\n"); printf("\n");
  }
  else if( (couplings == 2) && (dipoles == 1) )
  {
    dipolmom(nz, na, nb, nc, nd, qx, qy);
    readCouplings(nz, vab);
    printf("Couplings from File\n");
  }
  else if( (wscp == 1) && (couplings == 2) && (dipoles == 2) )  // Chl_a: 30°
  {
    char *in_dat_1 = "qy_30.dat";   

    readDipol2(nz, qy, in_dat_1);
    readCouplings(nz, vab);
  }
  else if( (wscp == 2) && (couplings == 2) && (dipoles == 2) )   // Chl_b: 36°  
  {
    char *in_dat_1 = "qy_36.dat";   // qy_30.dat ... qy_39.dat

    readDipol2(nz, qy, in_dat_1);
    readCouplings(nz, vab);
  }
  else
  {
    printf("Exit with couplings = %d, dipoles = %d\n", couplings, dipoles);
    exit(8);
  }
  if ( (E_BURN <= SPEKMIN) || (E_BURN >= SPEKMAX ) )
  {
    printf("\nFehler: Brennfrequenz außerhalb SPEKMIN und SPEKMAX\n");
    exit(8);
  }


//---------------------------------------------------------------------------------------------------------------


  calcJw(jw);                                       // Berechnung von J(w) eq 2.17,  [J(w)] = fs
  calcGt(gt, jw, in, out, plan_fwd);                // Berechnung von G(t) mit FFT,  [G(t)] = 
  calcCt(ct, jw, in, out, plan_fwd);                // Berechnung von C(t) mit FFT,  [C(t)] = fs-1
  calcCw(ct, cw, cw_re, cw_im, in, out, plan_bwd);  // Berechnung von C(w) mit FFT,  [C(w)] = cm-1


  lineshapeArray(gt, sideband, zeroline, in, out, plan_bwd);  // FastFourierTransf

  

  for(i=0; i<N; i++)  // BEGIN  Monte-Carlo-Loop 
  {
    calcLinearSpectrum( nz, i, anz, site, sigma, seed, vab, hamilt, eigval, eigvec, ev_h, qy, alpha, rotst, 
		        delta, dipst, rij, gam, tau, cw_re, cw_im, jw, gt, w_m0_sh, sideband, zeroline, 
		        spek_ab, sumspek_ab, spek_sb, arrayspek_sb, sumspek_sb, 
		        spek_zl, arrayspek_zl, sumspek_zl, spek_cd, sumspek_cd, in, out, plan_bwd );
       
    for(pigm_no=0; pigm_no<nz; pigm_no++) // Burn-Loop  (Pigment pigm_no wird gebrannt (Pigment-Nr, nicht Exciton-Nr!!))
    {
      calcHoleBurnSpectrum(nz, pigm_no, anz, wstep, site, sigma, seed, vab, vab_hb, hamilt, eigval, eigvec, eigval_hb, 
			   eigvec_hb, ev_h, qy, alpha, dipst, rij, gam, tau, cw_re, cw_im, jw, gt, w_m0_sh, w_hb_sh, 
			   sideband, zeroline, arrayspek_sb, arrayspek_zl, 
			   spek_ab, spek_sb, spek_zl, spek_hb, spek_hb_res, sumspek_hb, sumspek_hbab, in, out, plan_bwd); 
    
    } //End Burn-Loop
    
    

  } //END MC-Loop  


  //outputPigmExc(npig, nexc);
  

  spekDiff(anz, sumspek_hbab, sumspek_hb, sumspek_hbdif); 

  outputSpectrum(sumspek_ab,   "out_spec_absor.dat");  // absorbtion
  outputSpectrum(sumspek_sb,   "out_spec_sb_ab.dat");  // sideband
  outputSpectrum(sumspek_zl,   "out_spec_zl_ab.dat");  // zeroline
  outputSpectrum(sumspek_cd,   "out_spec_circd.dat");  // circular-dichrois
  outputSpectrum(sumspek_hbab, "out_spec_hburn.dat");  // pre-burn
  outputSpectrum(sumspek_hb,   "out_spec_hbabs.dat");  // post-burn
  outputSpectrum(sumspek_hbdif,"out_spec_hbdif.dat");  // hb-diff
  

	
  prgende=clock(); //CPU-Zeit am Ende des Programmes
  printf("\nLaufzeit %.2f Sekunden\n\n",(float)(prgende-prgstart) / CLOCKS_PER_SEC);

 

// ===== END Hauptprogramm ===================================================


// --- free allocated memory:


  free(num);  num  = NULL;

  free(x); x = NULL;
  free(y); y = NULL;
  free(z); z = NULL;

  free(dipst); dipst = NULL;
  free(site);  site = NULL;
  free(sigma); sigma = NULL;
  free(alpha); alpha = NULL;
  free(rotst); rotst = NULL;
  free(tau);   tau = NULL;
  free(w_m0_sh); w_m0_sh = NULL;
  free(w_hb_sh); w_hb_sh = NULL;

  free(jw); jw  = NULL;
 
  free(cw_im); cw_im = NULL;
  free(cw_re); cw_re = NULL;


  free(spek_ab);    spek_ab = NULL; 
  free(spek_sb);    spek_sb = NULL; 
  free(spek_zl);    spek_zl = NULL; 
  free(spek_cd);    spek_cd = NULL; 
  free(spek_hb);    spek_hb = NULL;
  free(spek_hb_res);spek_hb_res = NULL;
  free(sumspek_ab); sumspek_ab = NULL;
  free(sumspek_sb); sumspek_sb = NULL;
  free(sumspek_zl); sumspek_zl = NULL;
  free(sumspek_cd); sumspek_cd = NULL;
  free(sumspek_hb); sumspek_hb = NULL;
  free(sumspek_hbab); sumspek_hbab = NULL;
  free(sumspek_hbdif); sumspek_hbdif = NULL;


  // FFTW-complex-numbers

  fftw_free(ct); ct = NULL; 
  fftw_free(gt); gt = NULL;
  fftw_free(cw); cw = NULL;


// free vectors and matrices -------------------------------------------------------------


  nrow = nz;
  while(nrow--)
  {
    free(mg[nrow]);
    free(na[nrow]);
    free(nb[nrow]);
    free(nc[nrow]);
    free(nd[nrow]);
    free(qx[nrow]);
    free(qy[nrow]);
    free(vab[nrow]);
    free(vab_hb[nrow]);
    free(rij[nrow]);
    free(gam[nrow]);
    free(j_w[nrow]);
    free(n_w[nrow]);
    free(npig[nrow]);
    free(nexc[nrow]); 
    free(arrayspek_sb[nrow]); 
    free(arrayspek_zl[nrow]);

  }
  free(mg); mg = NULL;
  free(na); na = NULL;
  free(nb); nb = NULL;
  free(nc); nc = NULL;
  free(nd); nd = NULL;
  free(qx); qx = NULL;
  free(qy); qy = NULL;
  free(vab); vab = NULL;
  free(vab_hb); vab_hb = NULL;
  free(rij); rij = NULL;
  free(gam); gam = NULL;
  free(j_w); j_w = NULL;
  free(n_w); n_w = NULL;
  free(npig); npig = NULL;
  free(nexc); nexc = NULL;
  free(arrayspek_sb); arrayspek_sb = NULL;
  free(arrayspek_zl); arrayspek_zl = NULL;


  nrow = hn;
  while(nrow--)
  {
    free(delta[nrow]);
  }
  free(delta); delta = NULL;


  nrow = NGAM;
  while(nrow--)
  {
    ncol = NTAU;
    while(ncol--)
    {
      free(sideband[nrow][ncol]);
      free(zeroline[nrow][ncol]); 
    }
    free(sideband[nrow]);
    free(zeroline[nrow]);   
  }
  free(sideband);
  free(zeroline);   
  

  free(eigval); eigval = NULL;
  free(eigval_hb); eigval_hb = NULL;
  free(ev_h);   ev_h   = NULL;
  free(hamilt); hamilt = NULL;

  nrow = nz;
  while(nrow--)
  {
    free(eigvec[nrow]);
    free(eigvec_hb[nrow]);
  }
  free(eigvec); eigvec = NULL;
  free(eigvec_hb); eigvec_hb = NULL;



  // FFTW-plans

  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd); 
  fftw_free(in); 
  fftw_free(out);


  return (0);

}


//---------------------------------------------------------------------------------------------


    // Gauss-dressed-sticks der mittleren Site Energies:
    /*

   //double *ngau_abs = calloc( NGAU, sizeof(double));


    if(i==0)
    {
      printf("SiteEnergies: "); for(j=0;j<nz;j++) printf(" %6.0f", site[j]); printf("\n"); 
      printf("Eigenwerte:   "); for(j=0;j<nz;j++) printf(" %6.0f", nr_ew[j+1]); printf("\n");
      
      stick2gauss(nz, nr_ew, alpha, ngau_abs);
      normierung(NGAU, GAULENG, ngau_abs);
      outputGauss(nr_ew, ngau_abs, nz);
      outputStickSpek(nr_ew, alpha, nz);
    }
    */
