#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

#define EPS 1.0e-25

#define NMAX 5

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
  for (int i = 0; i < n; i++){
    y[i] = 0;
    for (int j = 0; j < n; j++){
      y[i] += a[i*n+j] * x[j];
    }
  }
}

double residual(int n, double *x, double *y){
	int i;
	double res = 0;

	for(i=0; i<n; i++){
		res += pow( (x[i]-y[i]), 2);
	}	
	return res;
}

double in(int n, double *x, double *y){
  double product = 0;
  for (int i = 0; i < n; i++){
    product += x[i] * y[i];
  }
  return product;
}


void cgbi_stab(int n, double *adim, double *bvec, double *xvec){
  int i, j, k, count = 0;
   double eps, oldeps, r, p, pp, rr;
   double alpha, beta, omega;
   double rap;

  double *pvec = new double [n]();
  double *rvec = new double [n]();
  double *r0vec = new double [n]();
  double *apvec = new double [n]();

  double *svec = new double [n]();

  double *asvec = new double [n]();


  for(i=0; i<n; i++){
    r0vec[i] = bvec[i];
    rvec[i] = pvec[i] = r0vec[i];
  }

  rr = in(n, r0vec, r0vec);

  if(rr == 0){
    exit(1);
  }

  while(1){
    oldeps = eps;
    eps = 0.0;

    //get apvec 
    map(n, adim, pvec, apvec);
    rr = in(n, r0vec, rvec);

    rap = in(n, r0vec, apvec);

    alpha = rr / rap;


    for(i=0; i<n; i++){
      svec[i] = rvec[i] - alpha * apvec[i];
    }


    //asvec = As
    map(n, adim, svec, asvec);

    omega = in(n, asvec, svec) / in(n, asvec, asvec);

    for(i=0; i<n; i++){
      xvec[i] += alpha * pvec[i] + omega * svec[i];
      rvec[i] = svec[i] - omega * asvec[i];
    }


    // double *bb = new double[n]();
    // map(n, adim, xvec, bb);

    //eps = residual(n, bvec, bb);

    for(i=0; i<n; i++){
      eps += fabs(rvec[i]);
    }
        
    //printf("%d %e   \n", count, eps);
    
    if(eps <= EPS){
      count++;
      break;
    }else if(isinf(eps)){
      break;
    }else if (isnan(eps)){
      break;
    }

    // for(i=0, rr=0, pp=0; i<n; i++){ 
    //   rr += rvec[i] * rvec[i];
    //   pp += pvec[i] * pvec[i];
    // }
    // beta = rr/pp;


    beta = (alpha / omega) * (in(n, r0vec, rvec) / rr);

    for(i=0; i<n; i++){
      pvec[i] = rvec[i] + beta * (pvec[i] - omega * apvec[i]);
    }
    count++;
  }

  printf("count: %d\n   eps: %e   oldeps: %e\n", count, eps, oldeps);
  free(pvec);free(rvec);free(apvec);
}

int main(int ac, char *av[]){

  int i,j,k;
  int n;
  double *adim;
  double *bvec;
  double *xvec;
  
  char data[30] = "./mt10";

  if(ac > 1){
    strcpy(data, av[1]);
  }

  printf("%s\n", data);

  readmats(data, n, adim, bvec, xvec);

  ////////////////////////////////////
  for(k=0; k<n; k++){
    for(j=0; j<n; j++){
      printf("%e	", adim[k*n+j]);
    }
    printf("\n");
  }
  printf("\n\nb\n");

  for(j=0; j<n; j++){
      printf("%lf ", bvec[j]);
  }
  printf("\n\nx\n");
  for(j=0; j<n; j++){
      printf("%lf ", xvec[j]);
  }
  printf("\n\n");
  ////////////////////////////////////

  cgbi_stab(n, adim, bvec, xvec);

  for(k=0; k<n; k++){
    printf("xvec[%2d] = %lf\n", k, xvec[k]);
  }

  free(adim);free(bvec);free(xvec);

  return 0;
}