/****************************************************************************

Julian Adolphs 2018

Das Program fitted die Site Energies mit einem Genetischen Algorithmus. 

Die Experimentellen Spektren werden eingelesen (Absorption, CD, LD) und mit einer Least-Square-Methode 
der Fitness-Wert des jeweiligen Spektrums berechnet. 
Es werden des Absorp, CD und LD Spektrum und die Ableitung des Absorp-Spektrums mit einbezogen.   

Die inhomogenen Spektren in dynamischer Theorie werden mit einer MonteCarlo-Methode berechnet.

Input:

 infile_koord.dat
 infile_koord_trim.dat
 infile_Exp_Spek.dat


Output:  

 fit_name.out       // Optimale Site Energies, in fit_fitnessFkt.c 

 out_spec_OD.dat    // Spektren für letzte Fit-Runde, wenn Konvergenz erreicht wurde, 
 out_spec_CD.dat    // entspricht das dem optimalen Satz von Site Energies 
 out_spec_LD.dat 
 out_spec_dxOD.dat

---

Files die zum Programm "spektrum" gehoeren:  spek_*.c spek_*.h

Files die zum Programm "fit" gehoeren:  fit_*.c fit_*.h

Fourier-Transformationen mit FFTW-Routinen 

Diagonalisierung der Hamilton-Matrix mit "Eispack": eispack.c eispack.h

Zufallszahlen: normal.h normal.c

****************************************************************************/
  

#include "spek_head.h"   // Makros & Definitionen
#include  "fit_head.h" 


int main(int argc, char *argv[ ]) 
{

  clock_t prgstart, prgende;  // Laufzeitmessung
  prgstart = clock();         // CPU-Zeit zu Beginn des Programmes

  int i=0, j=0, k=0, n=0, nz=0, trim_nz=0, hn=0, nrow=0, ncol=0; // m=0, pigm_no=0; 

  int bclnum = 8;        // 7 or 8 pigment variant of FMO

  int nmc_8 = N*BCL_8;   // Anteil MC-steps mit 8 pigmenten

  int h_var = 17;        // seed for random number generator
  int *seed  = &h_var;

  int species   = SPECIES;     // Chl_a or Chl_b-WSCP of Tep-FMO?!  

  char *in_dat = "infile_koord.dat";
  n = countInput(in_dat);

  char *in_dat_exp_1 = "infile_Exp_Spek.dat";  
  int nx = countInput(in_dat_exp_1);


  nz = n/NNATM;                  // Anzahl der chlorophylle, WSCP:2
  trim_nz = 3*nz;                // Anz Chls im Trimer
  hn = SQ(nz);
  printf("# Anzahl N-Atome:  %d ,  Anzahl Chl:  %d,  Temperatur:  %4.1lf K,  Winkel:  %4.1lf Grad\n", n, nz, TEMP, BETA);
  
  double wstep = TWOPI/(NT*DELTA*WZ2WFS);   // [cm-1]     
  int    anz   = (int) (SPEKDIFF/wstep + 3); 
  printf("# wstep:  %lf cm-1, anz:  %d\n", wstep, anz);

  if ( nx != anz )
  {
    printf("\nFehler: Experimentelle Spektren haben falsche Anzahl Stuetzstellen!\n");
    exit(8);
  }


  // geneticFit
  int nrekomb = ceil(NPOP*REKOMB_RATE);  // ceil = "Aufrunden"  
  int nrek_2  = ceil(NPOP*REK_2_RATE);
  int nmutat  = ceil(MUTAT_RATE*NPOP);   // ceil(x)  <=>  n => x, n int
  int nninsrt = ceil(NINSRT_RATE*NPOP);
  int n1 = nninsrt; 
  int n2 = nrekomb+nninsrt; 
  int n3 = nrek_2+nrekomb+nninsrt; 
  int n4 = nrekomb+nrek_2+nninsrt+nmutat;
  int npop = NPOP;




 // ------------- Declare & Allocate vectors  ---------------------------------

  int   *num = calloc( n, sizeof(int));
  int   *trim_num = calloc( trim_nz, sizeof(int) );

  double *x = calloc( n, sizeof(double));
  double *y = calloc( n, sizeof(double));
  double *z = calloc( n, sizeof(double));

  double *sym_ax =  calloc( 3, sizeof(double));    // Symmetrie-Achse fuer LD-Spek

  double *dipst   = calloc( nz, sizeof(double));
  double *site    = calloc( nz, sizeof(double));   // Site Energie Mittelwerte 
  double *sigma   = calloc( nz, sizeof(double));
  double *alpha   = calloc( nz, sizeof(double));   // alpha=|mu|^2 dipst^2 [in Debye^2] 
  double *rotst   = calloc( nz, sizeof(double));   // rotational strength for CD 
  double *ldcff   = calloc( nz, sizeof(double));   // LD-coefficient for LD 
  double *tau     = calloc( nz, sizeof(double));
  double *w_m0_sh = calloc( nz, sizeof(double));
  double *w_hb_sh = calloc( nz, sizeof(double));

  double *jw    = calloc( NT/2, sizeof(double));   // J(w)   jw[] = jw[0] ... jw[NT/2-1]

  double *cw_re = calloc( NT, sizeof(double));
  double *cw_im = calloc( NT, sizeof(double));

  double *spek_ab       = calloc( anz, sizeof(double)); // Homogenes Absorptions-Spektrum 
  double *spek_sb       = calloc( anz, sizeof(double)); // Homogenes Absorptions-Spektrum nur der SeitenBande
  double *spek_zl       = calloc( anz, sizeof(double)); // Homogenes Absorptions-Spektrum nur der NullNullLinie
  double *spek_cd       = calloc( anz, sizeof(double)); //
  double *spek_ld       = calloc( anz, sizeof(double)); // Homogenes LD-Spektrum 
  double *sumspek_ab    = calloc( anz, sizeof(double)); // Inhomogenes Absorptionsspektrum (=aufsummiert)
  double *sumspek_sb    = calloc( anz, sizeof(double)); // Inhomogenes AbsSpek der SB (=aufsummiert)
  double *sumspek_zl    = calloc( anz, sizeof(double)); // Inhomogenes AbsSpek der ZL (=aufsummiert)
  double *sumspek_cd    = calloc( anz, sizeof(double)); // 
  double *sumspek_ld    = calloc( anz, sizeof(double)); // Inhomogenes LD-spektrum (=aufsummiert)
  double *sumspek_dx    = calloc( anz, sizeof(double)); // Ableitung von sumspek_ab

  double *exp_ab       = calloc( nx, sizeof(double)); 
  double *exp_cd       = calloc( nx, sizeof(double)); 
  double *exp_ld       = calloc( nx, sizeof(double)); 
  double *exp_dx       = calloc( nx, sizeof(double));   // Ableitung OD'
  

  // genetic-fit
  double *fitness    = calloc( NPOP, sizeof(double));
  double *rank_fitn  = calloc( NPOP, sizeof(double));
  double *prob       = calloc( (NPOP+1), sizeof(double));
  double *vertlgsfkt = calloc( (NPOP+1), sizeof(double));
  double *nr_fitness = calloc( (NPOP+1), sizeof(double));

  int *nr_indx      = calloc( (NPOP+1), sizeof(int));
  int *nr_irank     = calloc( (NPOP+1), sizeof(int));
 
  int *irank      = calloc( NPOP, sizeof(int) );
  int *indx_sel   = calloc(1, sizeof(int));
  int *num_sel    = calloc(1, sizeof(int));
  int *indx_sel_1 = calloc(1, sizeof(int));
  int *num_sel_1  = calloc(1, sizeof(int));
  int *indx_sel_2 = calloc(1, sizeof(int));
  int *num_sel_2  = calloc(1, sizeof(int)); 
 



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

  // using FFTW-routines for fourier transforms
  plan_bwd = fftw_plan_dft_1d(NT, in, out, +1, FFTW_MEASURE);  // backward-transform
  plan_fwd = fftw_plan_dft_1d(NT, in, out, -1, FFTW_MEASURE);  // forward-transform


  //---- Matrizen -----------------------------------------------------------------------------------
  // Im FFTW-Manual wird von dieser Art Arrays zu konstruieren abgeraten, Schmaranz empfiehlt es so...

  double **parent    = calloc(npop, sizeof(double*)); //genetic-fit
  double **offspring = calloc(npop, sizeof(double*)); 
  nrow = npop;
  while(nrow--){   
    parent[nrow]    = calloc(nz, sizeof(double));
    offspring[nrow] = calloc(nz, sizeof(double));
  }

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

  double **trim_r = calloc( trim_nz, sizeof(double*));
  nrow = trim_nz;
  while(nrow--)
    trim_r[nrow] = calloc(3, sizeof(double));

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
    arrayspek_sb[nrow] = calloc( 4*anz+9, sizeof(double) ); 
    arrayspek_zl[nrow] = calloc( 4*anz+9, sizeof(double) );
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
  init(nz, site, sigma, dipst);                  
  koord(nz, x, y, z, na, nb, nc, nd);
  pigmCenters(nz, mg, nb, nd, delta, rij);   

  in_dat="infile_koord_trim.dat"; 
  readTRImerCoord(trim_num, trim_r, trim_nz, in_dat);
  trimSymAx(trim_nz/3, trim_r, sym_ax); 

  eingabeExpSpek(in_dat_exp_1, exp_ab, exp_cd, exp_ld, exp_dx, nx);   // Spektrum das gefittet werden soll einlesen

  maxnorm(nx, exp_ab);
  maxnorm(nx, exp_cd);
  maxnorm(nx, exp_ld);
  maxnorm(nx, exp_dx);
  //  outputSpectrum(exp_ab, anz, "out_spec_fit_OD.dat");
  //  outputSpectrum(exp_cd, anz, "out_spec_fit_CD.dat");
  //  outputSpectrum(exp_ld, anz, "out_spec_fit_LD.dat");
  //  outputSpectrum(exp_dx, anz, "out_spec_fit_dxOD.dat");



//---------------------------------------------------------------------------------------------------------------


  calcJw(jw);                                       // Berechnung von J(w) eq 2.17,  [J(w)] = fs
  calcGt(gt, jw, in, out, plan_fwd);                // Berechnung von G(t) mit FFT,  [G(t)] = 
  calcCt(ct, jw, in, out, plan_fwd);                // Berechnung von C(t) mit FFT,  [C(t)] = fs-1
  calcCw(ct, cw, cw_re, cw_im, in, out, plan_bwd);  // Berechnung von C(w) mit FFT,  [C(w)] = cm-1


  lineshapeArray(gt, sideband, zeroline, in, out, plan_bwd);  // FastFourierTransf

 
  for(k=0; k<NFIT; k++) // BEGIN k-Fit-Loop 
  {
    printf("Fit-Runde: %d \n", k);
 
    for(j=0; j<NPOP; j++) // BEGIN j-Loop ueber alle Chromosome 
    {  
      bclnum = 8; 
      parameterFit(site, nz, j, k, bclnum, dipst, parent);  // fit_parameter
      dipolmom(nz, na, nb, nc, nd, qx, qy);
      interaction(nz, qy, delta, vab, rij, dipst);
      
      for(i=0; i<nmc_8; i++)  // BEGIN  Monte-Carlo-Loop 
      {
	calcLinearSpectrum(nz, i, anz, site, sigma, seed, vab, hamilt, eigval, eigvec, ev_h, qy, alpha, rotst, ldcff,
			   delta, dipst, rij, gam, tau, cw_re, cw_im, jw, gt, w_m0_sh, sideband, zeroline, 
			   spek_ab, sumspek_ab, spek_sb, arrayspek_sb, sumspek_sb, 
			   spek_zl, arrayspek_zl, sumspek_zl, spek_cd, sumspek_cd, spek_ld, sumspek_ld, sym_ax, in, out, plan_bwd );
      } //END 8    
      
      bclnum = 7; 
      parameterFit(site, nz, j, k, bclnum, dipst, parent); 
      dipolmom(nz, na, nb, nc, nd, qx, qy);
      interaction(nz, qy, delta, vab, rij, dipst);

      for(i=nmc_8; i<N; i++)  // BEGIN 8
      {
	calcLinearSpectrum(nz, i, anz, site, sigma, seed, vab, hamilt, eigval, eigvec, ev_h, qy, alpha, rotst, ldcff,
			   delta, dipst, rij, gam, tau, cw_re, cw_im, jw, gt, w_m0_sh, sideband, zeroline, 
			   spek_ab, sumspek_ab, spek_sb, arrayspek_sb, sumspek_sb, 
			   spek_zl, arrayspek_zl, sumspek_zl, spek_cd, sumspek_cd, spek_ld, sumspek_ld, sym_ax, in, out, plan_bwd );
      } //END 7    
      

      derivative(anz, wstep, sumspek_ab, sumspek_dx); 
      
      maxnorm(anz, sumspek_ab);
      maxnorm(anz, sumspek_cd);
      maxnorm(anz, sumspek_ld);
      maxnorm(anz, sumspek_dx);

      if( k==(NFIT-1) && ( j==0 ) )  // plot spektrum for best-fitness Chromosome (usually no 1) 
      {
	 outputSpectrum(sumspek_ab, anz, "out_spec_OD.dat");  
	 outputSpectrum(sumspek_cd, anz, "out_spec_CD.dat");  
	 outputSpectrum(sumspek_ld, anz, "out_spec_LD.dat"); 
	 outputSpectrum(sumspek_dx, anz, "out_spec_dxOD.dat");  
      }
      
      // edit, which part of spectrum ist fitted
      diffExpSim(j, anz, wstep, fitness, nr_fitness, sumspek_ab, sumspek_cd, sumspek_ld, sumspek_dx, exp_ab, exp_cd, exp_ld, exp_dx); 

      reset(sumspek_ab, anz); 
      reset(sumspek_cd, anz); 
      reset(sumspek_ld, anz);
      reset(sumspek_dx, anz);

    } //END j-Loop ueber Chromosomen  

    
    indexxing(npop, nr_fitness, nr_indx);
    ranking(npop, nr_indx, nr_irank);      
    sorting(npop, nr_fitness); 

    verteilungsfkt(fitness, nr_fitness, rank_fitn, vertlgsfkt, prob, irank, nr_irank, parent);  

    geneticOperate(n1, n2, n3, n4, nz, i, k, irank, num_sel, num_sel_1, num_sel_2, 
    		   indx_sel, indx_sel_1, indx_sel_2, fitness, vertlgsfkt, parent, offspring);
   

  } //END k-Loop ueber Fit-steps
 
	
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
  free(ldcff); ldcff = NULL;
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
  free(sumspek_ab); sumspek_ab = NULL;
  free(sumspek_sb); sumspek_sb = NULL;
  free(sumspek_zl); sumspek_zl = NULL;
  free(sumspek_cd); sumspek_cd = NULL;
  free(sumspek_dx); sumspek_dx = NULL;

  free(exp_ab); exp_ab = NULL; 
  free(exp_ld); exp_ld = NULL; 
  free(exp_cd); exp_cd = NULL; 
  free(exp_dx); exp_dx = NULL; 

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



	/* HB-Spek unguenstig fuer Fit wegen der Normierungsprobleme

	for(pigm_no=0; pigm_no<nz; pigm_no++) // Burn-Loop  (Pigment pigm_no wird gebrannt (Pigment-Nr, nicht Exciton-Nr!!))
	{
	  calcHoleBurnSpectrum(nz, pigm_no, anz, wstep, site, sigma, seed, vab, vab_hb, hamilt, eigval, eigvec, eigval_hb, 
			   eigvec_hb, ev_h, qy, alpha, dipst, rij, gam, tau, cw_re, cw_im, jw, gt, w_m0_sh, w_hb_sh, 
			   sideband, zeroline, arrayspek_sb, arrayspek_zl, 
			   spek_ab, spek_sb, spek_zl, spek_hb, spek_hb_res, sumspek_hb, sumspek_hbab, in, out, plan_bwd); 
	} //End Burn-Loop
	*/

      //spekDiffHB(anz, sumspek_hbab, sumspek_hb, sumspek_hbdif); // calculation of HB-Spektrum


/*
  char *in_dat_exp_2 = "exp_spek_absor.dat"; //"exp_spek_hburn.dat";  //"test_spek.dat";            
  int nx_2 = countInput(in_dat_exp_2);

  if ( nx != nx_2 )
  {
    printf("\nFehler: Experimentelle Spektren haben unterschiedlich viele Stuetzstellen!\n");
    exit(8);
  }
*/


 /*  
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
  */




/*
  if ( (E_BURN <= SPEKMIN) || (E_BURN >= SPEKMAX ) )
  {
    printf("\nFehler: Brennfrequenz außerhalb SPEKMIN und SPEKMAX\n");
    exit(8);
  }
*/




  /*
  outputSpectrum(sumspek_ab,   "out_spec_absor.dat");  // absorption
  outputSpectrum(sumspek_sb,   "out_spec_sb_ab.dat");  // sideband
  outputSpectrum(sumspek_zl,   "out_spec_zl_ab.dat");  // zeroline
  outputSpectrum(sumspek_cd,   "out_spec_circd.dat");  // circular-dichrois
  outputSpectrum(sumspek_hbab, "out_spec_hburn.dat");  // pre-burn
  outputSpectrum(sumspek_hb,   "out_spec_hbabs.dat");  // post-burn
  outputSpectrum(sumspek_hbdif,"out_spec_hbdif.dat");  // hb-diff
  */

