#include "invert.h"

int main(int argc, char **argv)
{
  CPREC **a,**acopy,**b;
  int n = 7,i,j;
  FILE *in;
  if(argc>1)
    sscanf(argv[1],"%i",&n);
  a = nr_matrix(1,n,1,n);
  acopy = nr_matrix(1,n,1,n);
  b = nr_matrix(1,n,1,1);
  
  in = fopen("a.dat","r");
  for(i=1;i<=n;i++){
    b[i][1] = 1.0;
    for(j=1;j<=n;j++){
      fscanf(in,FMT_1,&(acopy[i][j]));
    }
    
  }
  fclose(in);

  
  for(i=1;i<=n;i++){
    b[i][1] = 1.0;
    for(j=1;j<=n;j++)
      a[i][j]=acopy[i][j];
  }
  gaussj(a,n,b,1);



  for(i=1;i<=n;i++){
    b[i][1] = 1.0;
    for(j=1;j<=n;j++)
      a[i][j]=acopy[i][j];
  }
  lu_solve(a,n,b,1);

  return 0;

}
