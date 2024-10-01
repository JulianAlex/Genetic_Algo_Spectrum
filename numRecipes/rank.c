/* ju-changed from float 2 double */

void rank(unsigned long n, unsigned long indx[], unsigned long irank[])

     /* Given "indx[1..n]" as output from the routine "indexx", returns an array
        "irank[1..n]", the corresponding table of ranks */

{
	unsigned long j;

	for (j=1;j<=n;j++) irank[indx[j]]=j;
}
/* (C) Copr. 1986-92 Numerical Recipes Software '>'!^,. */
