#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

void readmats(char *fnm, int &n, double *&a, double *&b, double *&x){
  std::ifstream ifs;
  ifs.open(fnm);
  ifs >> n;
  a = new double [n*n];
  b = new double [n];
  x = new double [n]();

  
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      ifs >> a[i*n+j];
    }
    ifs >> b[i];
  }
  ifs.close();
}

int main(int ac, char *av[]){

  int i,j,k;
  int n;
  double *adim;
  double *bvec;
  double *xvec;
  
  char data[30] = "./mt1000";

  printf("Data set: %s\n", data);

  //read File
  readmats(data, n, adim, bvec, xvec);


 puts("A");
  for(k=0; k<n; k++){
    for(i=0; i<n; i++){
      printf("%lf  ", adim[k*n+i]);
    }
    printf("\n");
  }
  printf("\nb\n");
  for(i=0; i<n; i++){
      printf("%lf  ", bvec[i]);
    }
  puts("");


  free(adim);free(bvec);free(xvec);

  return 0;
}
