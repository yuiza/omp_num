#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

#define EPS 1.0e-15

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

void map(int n, double *a, double *x, double *b){
  for (int i = 0; i < n; i++){
    b[i] = 0;
    for (int j = 0; j < n; j++){
      b[i] += a[i*n+j] * x[j];
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

void gs(int n, double *adim, double *bvec, double *xvec){
  int i, j, k, count = 1;
  double eps;
  double old;
  double *bdash = new double[n]();

  while(1){
    eps = 0.0;
  
    //map(n, adim, bvec, xvec);
	
    for(k=0; k<n; k++){
 	    old = xvec[k];    
      xvec[k] = bvec[k];
      for(i=0; i<n; i++){
        if(i != k)
          xvec[k] -= adim[k*n+i] * xvec[i];
      }
      xvec[k] =  xvec[k] / adim[k*n+k];
    }

    map(n, adim, xvec, bdash);


    eps = residual(n, bvec, bdash);
    
    printf("%d %e   \n", count, eps);

    if(eps < EPS){
      count++;
      break;
    // }else if(isinf(eps)){
    //   puts("Inf");
    //   break;
    }else if(isnan(eps)){
      puts("Nan");
      break;
    }
    count++;
  }
}

void gs2(int n, double *adim, double *bvec, double *xvec){
   int i, j, k, count = 1;
   double eps, tmp;

   //gauss-s
   while(1){
      eps = 0.0;
      for(k=0; k<n; k++){
         // count++;
         tmp = xvec[k];
         xvec[k] = bvec[k];
         for(i=0; i<n; i++){
            if(i != k)
               xvec[k] -= adim[k*n+i] * xvec[i];
         }   
         xvec[k] =  xvec[k] / adim[k*n+k];   

         //εチェック
         eps += fabs(tmp - xvec[k]);
      }
      printf("%d %e   \n", count, eps);
      if(eps < EPS){
         count++;
         // printf("\nEND  %e ", eps);
         break;
      }else if(isinf(eps)){
         printf("\nERROR : eps inf\n");
         break;
      }else if(isnan(eps)){
         printf("\nERROR : eps nan\n");
         break;
      }
      count++;
   }
   //printf("%d\n", count);
}

double residua1(int n, double *r){
        int i;
        double res;

        for(i=0; i<n; i++){
                res += fabs(r[i]);
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


void cg(int n, double *adim, double *bvec, double *xvec){
  int i, j, k, count = 0;
   double eps, oldeps, alpha, beta, r, p, pp, rr;

  double *pvec = new double [n]();
  double *rvec = new double [n]();
  double *apvec = new double [n]();

  for(i=0; i<n; i++){
    rvec[i] = bvec[i];
    pvec[i] = rvec[i];
  }

  while(1){
    oldeps = eps;
    eps = 0.0;

    for(i=0; i<n; i++){
      apvec[i] = 0;
      for(j=0; j<n; j++){
        apvec[i] += adim[i*n+j] * pvec[j];
      }
    }

    for(i=0, r=0, p=0; i<n; i++){
      r += rvec[i] * rvec[i];
      p += pvec[i] * apvec[i];
    }
    alpha = r / p;

    for(i=0; i<n; i++){
      xvec[i] += alpha * pvec[i];
      rvec[i] -= alpha * apvec[i];
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

    for(i=0, rr=0, pp=0; i<n; i++){ 
      rr += rvec[i] * rvec[i];
      pp += pvec[i] * pvec[i];
    }
    beta = rr/pp;

    for(i=0; i<n; i++){
      pvec[i] = rvec[i] + beta * pvec[i];
    }
    count++;
  }
  // while(1){
  //   eps = 0.0;

  //   map(n, adim, pvec, apvec);
  //   double alpha = in(n, rvec,rvec) / in(n, pvec, apvec);

  //   for(i=0; i<n; i++){
  //     xvec[i] += alpha * pvec[i];

  //     rvec[i] -= alpha * apvec[i];
  //   }

  //   eps = residua1(n, rvec); 
    
  //   oldeps = eps;

  //   //printf("eps: %e\n", oldeps);

  //   if(eps <= EPS){
  //     count++;
  //     break;
  //   }else if(isinf(eps)){
  //     printf("eps: %e", oldeps);
  //     break;
  //   }else if (isnan(eps)){
  //     printf("eps: %e", oldeps);
  //     break;
  //   }

  //   double beta = in(n, rvec, rvec) / in(n, pvec, pvec);

  //   for(i=0; i<n; i++){
  //     pvec[i] = rvec[i] + beta * pvec[i];
  //   }
  //   count++;
  // }
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

  gs2(n, adim, bvec, xvec);
  //cg(n, adim, bvec, xvec);

  for(k=0; k<n; k++){
    printf("xvec[%2d] = %lf\n", k, xvec[k]);
  }

  free(adim);free(bvec);free(xvec);

  return 0;
}