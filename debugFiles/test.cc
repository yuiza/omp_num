#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

#include "timer.h"

#define EPS 1.0e-25
#define NMAX 10000

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

 double *tmp = new double [n]();

#pragma omp parallel for 
 for (int i = 0; i < n; i++){
#pragma omp parallel for
    for (int j = 0; j < n; j++){
      tmp[i] += a[i*n+j] * x[j];
    }
    y[i] = tmp[i];
  }
  free(tmp);
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
  int count;
  double eps, rr, rap;
  double alpha, beta, omega;
  double cnt;
  double ass, asas, r0r;

  double *pvec = new double [n]();
  double *rvec = new double [n]();
  double *r0vec = new double [n]();
  double *apvec = new double [n]();
  double *svec = new double [n]();
  double *asvec = new double [n]();

  double ts = dgettime();

#pragma omp parallel
{//parallel start

#pragma omp for
  for(int i=0; i<n; i++){
    r0vec[i] = bvec[i];
    rvec[i] = pvec[i] = r0vec[i];
  }
  count = 0;
  while(count < NMAX){
    eps = 0.0;

    //get apvec 
    //map(n, adim, pvec, apvec);
    //double *tmp = new double [n]();
#pragma omp for
 for (int i = 0; i < n; i++){
    apvec[i] = 0;
    for (int j = 0; j < n; j++){
      apvec[i] += adim[i*n+j] * pvec[j];
    }
    //apvec[i] = tmp[i];
  }


    //(r0, r)
    //rr = in_product(n, r0vec, rvec);
    rr = 0;
#pragma omp for reduction(+:rr)
  for (int i = 0; i < n; i++){
    rr += r0vec[i] * rvec[i];
  }


    //(r0, ap)
    //rap = in_product(n, r0vec, apvec);
    rap = 0;
#pragma omp for reduction(+:rap)
  for (int i = 0; i < n; i++){
    rap += r0vec[i] * apvec[i];
  }


    alpha = rr / rap;

    //svec
#pragma omp for 
    for(int i=0; i<n; i++){
      svec[i] = rvec[i] - alpha * apvec[i];
    }

    //asvec = As
    //map(n, adim, svec, asvec);

#pragma omp for
 for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      asvec[i] += adim[i*n+j] * svec[j];
    }
    //y[i] = tmp[i];
  }
  
  ass = 0;
#pragma omp for reduction(+:ass)
  for (int i = 0; i < n; i++){
    ass += asvec[i] * svec[i];
  }
  asas = 0;
#pragma omp for reduction(+:asas)
  for (int i = 0; i < n; i++){
    asas += asvec[i] * asvec[i];
  }

    //omega = in_product(n, asvec, svec) / in_product(n, asvec, asvec);
  omega = ass / asas;

#pragma omp for
    for(int i=0; i<n; i++){
      xvec[i] += alpha * pvec[i] + omega * svec[i];
      rvec[i] = svec[i] - omega * asvec[i];
    }

//#pragma omp for reduction(+:eps)
    for(int i=0; i<n; i++){
      eps += fabs(rvec[i]);            
    }

    if(eps < EPS){
      count++;
      printf("count: %d\neps: %e\n",cnt, eps);
      count = NMAX+1;
      continue;
      //count++;
      //break;
    }
    // else if(isinf(eps)){
    //   break;
    // }else if (isnan(eps)){
    //   break;
    // }

   r0r = 0;
#pragma omp for reduction(+:r0r)
  for (int i = 0; i < n; i++){
    r0r += r0vec[i] * rvec[i];
  }

  beta = (alpha / omega) * (r0r / rr);

    //beta = (alpha / omega) * (in_product(n, r0vec, rvec) / rr);

#pragma omp for
    for(int i=0; i<n; i++){
      pvec[i] = rvec[i] + beta * (pvec[i] - omega * apvec[i]);
    }


    count++;
  }
  
}//parallel end
  double te = dgettime();

  printf("TIME: %f\n", te-ts);

  //ループ数，残差
  //printf("count: %d\neps: %e\n", count, eps);

  //free忘れずに
  free(pvec);free(rvec);free(apvec);free(r0vec);free(svec);free(asvec);
}

double residual(int n, double *x, double *y){
  double res = 0;

  for(int i=0; i<n; i++){
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
puts("Read File Success");
  //BiCG Stab
  cgbi_stab(n, adim, bvec, xvec);

  //check x answer
  //for(int k=0; k<n; k++){
  //  printf("xvec[%2d] = %lf\n", k, xvec[k]);
  //}

  //check b answer
  check_answer(n, adim, xvec, bvec);


  free(adim);free(bvec);free(xvec);

  return 0;
}
