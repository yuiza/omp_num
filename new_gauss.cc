#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

#include "timer.h"

#define EPS 1.0e-15
#define NMAX 1000000

//ファイル読み込み
void readmats(char *fnm, int &n, double *&a, double *&b, double *&x){
  std::ifstream ifs;
  ifs.open(fnm);
  ifs >> n;
  a = new double [n*n];
  b = new double [n];
  x = new double [n];

  memset(x, 0, sizeof(double)*n);
	
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
		  ifs >> a[i*n+j];
		}
		ifs >> b[i];
  }
  ifs.close();
}

void gs(int n, double *adim, double *bvec, double *xvec){
        double *acp = new double [n*n];// a copy
        double *bcp = new double [n]; //b copy
	//double /*alpha = 0, */sum = 0;
        //int i, k, j;


	//memcpy(adim, a, sizeof(a));
#pragma omp parallel for 
	for(int i=0; i<n; i++){
	  for(int j=0; j<n; j++){
		acp[i+n*j] = adim[i+n*j];
	  }
	  bcp[i] = bvec[i];
	}
////////////////////////////////////////////////////////////////	
	double ts = dgettime();

        for(int k=0; k<n-1; k++){
#pragma omp parallel for firstprivate(k)
                for(int i=k+1; i<n; i++){
			double alpha =  -acp[i*n+k] / acp[k*n+k];
                        for(int j=0; j<n; j++){
				acp[i*n+j] += alpha * acp[k*n+j];
                        }
                        bcp[i] += alpha*bcp[k];
                }
        }

        for(int k=n-1; k>=0; k--){
		double sum = bcp[k];
#pragma omp parallel for firstprivate(k) reduction(-:sum)
                for(int j=k+1; j<n; ++j){
			sum -= acp[k*n+j] * xvec[j];
                }
                xvec[k] = sum/acp[k*n+k];
        }

	double te = dgettime();

	printf("TIME: %lf\n", te-ts);

	free(acp);
}

//残差
double residual(int n, double *x, double *y){
  int i;
  double res = 0;

  for(i=0; i<n; i++){
    res += pow( (x[i]-y[i]), 2);
  } 
  return res;
}

//検算部
void check_answer(int n, double *a, double *x, double *b){
  double *bb = new double [n]();
  double res;

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      bb[i] += a[i*n+j] * x[j];
    }
  }

  res = residual(n, b, bb);

  printf("res %e\n", res);

  free(bb);
}

int main(int ac, char *av[]){

  int i,j,k;
  int n;
  double *adim;
  double *bvec;
  double *xvec;
  
  char data[30] = "/share/data/matrix/mt1000";

  if(ac > 1){
    strcpy(data, av[1]);
  }
  printf("%s\n", data);

  //read File
  readmats(data, n, adim, bvec, xvec);

  puts("read Success");

  //BiCG Stab
  //cgbi_stab(n, adim, bvec, xvec);

  gs(n, adim, bvec, xvec);

  //check x answer
  //for(k=0; k<n; k++){
  //  printf("xvec[%2d] = %lf\n", k, xvec[k]);
  //}

  //check b answer
  check_answer(n, adim, xvec, bvec);


  free(adim);free(bvec);free(xvec);

  return 0;
}
