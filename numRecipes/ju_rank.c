void ju_rank(int n, int *indx, int *irank)

     /* Given "indx[1..n]" as output from the routine "indexx", returns an array
        "irank[1..n]", the corresponding table of ranks */

{
	int j;

	for (j=1;j<=n;j++) irank[indx[j]]=j;
}
/* (C) Copr. 1986-92 Numerical Recipes Software '>'!^,. */
