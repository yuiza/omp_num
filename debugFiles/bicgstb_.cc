#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

#define EPS 1.0e-2
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

//行列積
void map(int n, double *a, double *x, double *y){
#pragma omp parallel for
  for (int i = 0; i < n; i++){
    y[i] = 0;
    for (int j = 0; j < n; j++){
      y[i] += a[i*n+j] * x[j];
    }
    //printf("%f\n", y[i]);
  }
}

//内積
double in_product(int n, double *x, double *y){
  double product = 0;
#pragma omp parallel for reduction(+:product)
  for (int i = 0; i < n; i++){
    product += x[i] * y[i];
  }
  return product;
}

//norm
double norm(int n, double *vec){
	double norm = 0;
#pragma omp parallel for reduction(+:norm)
	for(int i=0; i<n; i++){
		norm += fabs(vec[i]);
	}
	return norm;
}

void cgbi_stab(int n, double *adim, double *bvec, double *xvec){
  int i, count;
  double eps, rr, rap, rrk;
  double alpha, beta, omega;

  double *axvec = new double [n]();
  double *pvec  = new double [n]();
  double *rvec  = new double [n]();
  double *r0vec = new double [n]();
  double *apvec = new double [n]();
  double *svec  = new double [n]();
  double *asvec = new double [n]();

  map(n, adim, xvec, axvec);

#pragma omp parallel for
  for(i=0; i<n; i++){
    r0vec[i] = bvec[i] - axvec[i];
    rvec[i]  = r0vec[i];
    pvec[i]  = rvec[i];
  }

  rr = in_product(n, r0vec, rvec);

  if(rr == 0){
    printf("rr == 0\n");
    exit(1);
  }

  count = 0;
  while(count < NMAX){
    eps = 0.0;

    //get apvec 
    map(n, adim, pvec, apvec);

    //(r0, ap)
    rap = in_product(n, r0vec, apvec);

    alpha = rr / rap;

    //svec
#pragma omp parallel for
    for(i=0; i<n; i++){
      svec[i] = rvec[i] - alpha * apvec[i];
    }

    //asvec = A * svec
    map(n, adim, svec, asvec);

    omega = in_product(n, svec, asvec) / in_product(n, asvec, asvec);

#pragma omp parallel for
    for(i=0; i<n; i++){
      xvec[i] = xvec[i] + alpha * pvec[i] + omega * svec[i];
      rvec[i] = svec[i] - omega * asvec[i];
    }

	eps = norm(n, rvec) / norm(n, bvec);

    if(eps <= EPS){
      count++;
      printf("count: %d\neps: %e\n", count, eps);
      count = NMAX +1;
      continue;
    }
     else if(isinf(eps)){
      count++;
      printf("INF\ncount: %d\neps: %e\n", count, eps);
      count = NMAX +1;
      continue;
     }else if (isnan(eps)){
      count++;
      printf("NAN\ncount: %d\neps: %e\n", count, eps);
      count = NMAX +1;
      continue;
     }

    //beta = alpha / omega * rrk+1 / rr
    rrk = in_product(n, r0vec, rvec);

    beta = (alpha / omega) * (rrk / rr);

    rr = rrk;

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

// 残差res
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

// 対象行列生成　51..151..15の行列    
void init(int n, double *&a, double *&b, double *&x){
  int i,j;

  a = new double [n*n]();
  b = new double [n]();
  x = new double [n]();

  for(i=0; i<n; i++){
      for(j=0; j<n; j++){
        if(i==0){
          a[i*n+0] = 5;
          a[i*n+1] = 1;
        }else if(i == n-1){
          a[i*n+i] = 5;
          a[i*n+i-1] = 1;
        }else if(i == j){
          a[i*n+j] = 5;
          a[i*n+j-1] = 1;
          a[i*n+j+1] = 1;
        }
      }
      if(i == 0 || i == n-1) b[i] = 6;
      else                   b[i] = 7;
  }
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

//対象行列用
  // n = 10;
  // if(ac > 1){
  //     n = atoi(av[1]);
  // }
  // printf("%d\n", n);
//ここまで

  //read File
  readmats(data, n, adim, bvec, xvec);

  // 対象行列生成用
  // init(n, adim, bvec, xvec);

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
