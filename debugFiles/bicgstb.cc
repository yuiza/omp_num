#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

#define EPS 1.0e-25
#define NMAX 1000000

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

void map(int n, double *a, double *x, double *y){
#pragma omp parallel for
  for (int i = 0; i < n; i++){
    y[i] = 0;
    for (int j = 0; j < n; j++){
      y[i] += a[i*n+j] * x[j];
    }
  }
}

double in_product(int n, double *x, double *y){
  double product = 0;
#pragma omp parallel for reduction(+:product)
  for (int i = 0; i < n; i++){
    product += x[i] * y[i];
  }
  return product;
}

//参考サイト：http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95
void cgbi_stab(int n, double *adim, double *bvec, double *xvec){
  int i, j, k, count;
  double eps, rr, rap;
  double alpha, beta, omega;

  double *pvec = new double [n]();
  double *rvec = new double [n]();
  double *r0vec = new double [n]();
  double *apvec = new double [n]();
  double *svec = new double [n]();
  double *asvec = new double [n]();

#pragma omp parallel for
  for(i=0; i<n; i++){
    r0vec[i] = bvec[i];
    rvec[i] = r0vec[i];
    pvec[i] = rvec[i];
  }

  count = 0;
  while(count < NMAX){
    eps = 0.0;

    //get apvec 
    map(n, adim, pvec, apvec);

    //(r0, r)
    rr = in_product(n, r0vec, rvec);

    //(r0, ap)
    rap = in_product(n, r0vec, apvec);

    alpha = rr / rap;

    //svec
#pragma omp parallel for
    for(i=0; i<n; i++){
      svec[i] = rvec[i] - alpha * apvec[i];
    }

    //asvec = As
    map(n, adim, svec, asvec);

    omega = in_product(n, asvec, svec) / in_product(n, asvec, asvec);

#pragma omp parallel for
    for(i=0; i<n; i++){
      xvec[i] += alpha * pvec[i] + omega * svec[i];
      rvec[i] = svec[i] - omega * asvec[i];
    }


#pragma omp parallel for reduction(+:eps)
    for(i=0; i<n; i++){
      eps += fabs(rvec[i]);            
    }

    if(eps < EPS){
      count++;
      printf("count: %d\neps: %e\n", count, eps);
      count = NMAX +1;
      continue;
      //break;
    }
    // else if(isinf(eps)){
    //   break;
    // }else if (isnan(eps)){
    //   break;
    // }

    beta = (alpha / omega) * (in_product(n, r0vec, rvec) / rr);

#pragma omp parallel for
    for(i=0; i<n; i++){
      pvec[i] = rvec[i] + beta * (pvec[i] - omega * apvec[i]);
    }

    count++;
  }

  //printf("count: %d\neps: %e\n", count, eps);

  //free忘れずに
  free(pvec);free(rvec);free(apvec);free(r0vec);free(svec);free(asvec);
}

double residual(int n, double *x, double *y){
  int i;
  double res = 0;

  for(i=0; i<n; i++){
    res += pow( (x[i]-y[i]), 2);
  } 
  return res;
}

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
  
  char data[30] = "./mt1000";

  if(ac > 1){
    strcpy(data, av[1]);
  }
  printf("%s\n", data);


  //read File
  readmats(data, n, adim, bvec, xvec);

  puts("read Success");

  //BiCG Stab
  cgbi_stab(n, adim, bvec, xvec);

  //check x answer
  for(k=0; k<n; k++){
    printf("xvec[%2d] = %lf\n", k, xvec[k]);
  }

  //check b answer
  check_answer(n, adim, xvec, bvec);


  free(adim);free(bvec);free(xvec);

  return 0;
}
