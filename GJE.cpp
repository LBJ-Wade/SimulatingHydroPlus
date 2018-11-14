#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

int gaussj(double **a, int n,double **b, int m)
//adapted from NR
//a[1...n][1...n] input matrix, b[1..n][1..m] r.h.s. vector
{

  //int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,temp;

  //test data
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      if (isnan(a[i][j])!=0)
	return 1;

  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      if (isnan(b[i][j])!=0)
	return 1;

  
  
  //indxc=VEC_I_INT(1,n);
  //indxr=ivector(1,n);
  //ipiv=ivector(1,n);
  int indxc[n+1];
  int indxr[n+1];
  int ipiv[n+1];


  for (j=0;j<n;j++) ipiv[j]=0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if (ipiv[j] !=1)
	for (k=0;k<n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k])>= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  }
	}
    ++(ipiv[icol]);
    if (irow!=icol) {
      for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
      for (l=0;l<m;l++) SWAP (b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) printf("gaussj: singular matrix\n");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    for (l=0;l<m;l++) b[icol][l] *= pivinv;
    
    for (ll=0;ll<n;ll++)
      if (ll != icol)
	{
	  dum=a[ll][icol];
	  a[ll][icol]=0.0;
	  for (l=0;l<n;l++) a[ll][l]-=a[icol][l]*dum;
	  for (l=0;l<m;l++) b[ll][l]-=b[icol][l]*dum;
	}
  }

  for (l=n-1;l>=0;l--) {
    if (indxr[l] !=indxc[l])
      for (k=0;k<n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
 
  return 0;
}

/*

double **a;
a = new double*[2];
for (int i=0;i<2;i++) a[i] = new double[2];
a[0][0]=1.0;
a[0][1]=3.0;
a[1][0]=2.0;
a[1][1]=1.0;

double **b;
b = new double*[2];
for (int i=0;i<2;i++) b[i] = new double[1];
b[0][0]=5.0;
b[1][0]=3.0;

gaussj(a,2,b,1);

cout << "b[0][0]: " << b[0][0];
cout << "b[1][0]: " << b[1][0] << endl;
*/
