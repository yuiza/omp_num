#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <fstream>
#include <omp.h>

#define EPS 1.0e-28

#define NMAX 5

void readmats(char *fnm, int &n, double *&a, double *&b){
	std::ifstream ifs;
	ifs.open(fnm);
	ifs >> n;
	a = new double [n*n];
	b = new double [n];
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			ifs >> a[i*n+j];
		}
		ifs >> b[i];
	}
	ifs.close();
}

void map(int size, double *a, double *x, double *y){
        for (int i = 0; i < size; i++){
                y[i] = 0;
                for (int j = 0; j < size; j++){
                        y[i] += a[i*size+j] * x[j];
                }
        }
}

double residual(int n, double *old, double *vec){
	int i;
	double res = 0;

	for(i=0; i<n; i++){
		res += pow( (old[i]-vec[i]) , 2);
	}	
	return res;
}

void gs(int n, double *adim, double *bvec, double *xvec){
   int i, j, k, count = 1;
   double eps;
   double *old = new double [n]();

   while(1){
      eps = 0.0;
      
      map(n, adim, bvec, xvec);
	
      for(k=0; k<n; k++){
         old[k] = xvec[k];
	//get new xvec 
	//map(n, adim, bvec, xvec);
         //xvec[k] = bvec[k];
         for(i=0; i<n; i++){
            if(i != k)
//               xvec[k] -= adim[k][i] * xvec[i];
               xvec[k] -= adim[k*n+i] * xvec[i];
         }
         xvec[k] =  xvec[k] / adim[k*n+k];

//        xvec * adim = res
//        eps += pow(old - new)

//        double eps = residual(n, old, xvec);

//        eps += fabs(old - xvec[k]);

     }
     map(n, adim, bvec, xvec);
     printf("%d %e   \n", count, eps);

      eps = residual(n, old, xvec);

      if(eps < EPS){
         count++;
         break;
      }else if(isinf(eps)){
        puts("Inf");
        break;
      }else if(isnan(eps)){
        puts("Nan");
        break;
      }
      count++;
   }
   printf("count: %d\n", count);
}

double residua1(int n, double *r){
        int i;
        double res;

        for(i=0; i<n; i++){
                res += fabs(r[i]);
        }
        return res;
}

void map2(int size, double *a, double *x, double *y){
	for (int i = 0; i < size; i++){
		y[i] = 0;
		for (int j = 0; j < size; j++){
			y[i] += a[i*size+j] * x[j];
		}
	}
}

double in(int size, double *x, double *y){
	double product = 0;
	for (int i = 0; i < size; i++){
		product += x[i] * y[i];
	}
	return product;
}


void cg(int n, double *adim, double *bvec, double *xvec){
	int i, j, k, count = 1;
   	double eps, alpha, beta, r, p, pp, rr, oldeps;

	double *pvec = new double [n]();
	double *rvec = new double [n]();
	double *apvec = new double [n]();

//	apvec = (double *)malloc(sizeof(double)*n);
//	pvec  = (double *)malloc(sizeof(double)*n);
//	rvec  = (double *)malloc(sizeof(double)*n);

	for(i=0; i<n; i++){
		rvec[i] = bvec[i];
		pvec[i] = rvec[i];
	}

	while(1){
		eps = 0.0;
		//for(i=0; i<n; i++){
		//	apvec[i] = 0;
		//	for(j=0; j<n; j++){
		//		apvec[i] += adim[i*n+j] * pvec[j];
		//	}
		//}
		map2(n, adim, pvec, apvec);
		//for(i=0, r=0, p=0; i<n; i++){
		//	r += rvec[i] * rvec[i];
		//	p += pvec[i] * apvec[i];
		//}
		//alpha = r / p;
		alpha = in(n, rvec,rvec) / in(n, pvec, apvec);

		for(i=0; i<n; i++){
			xvec[i] += alpha * pvec[i];
			rvec[i] -= alpha * apvec[i];
		}

      		double eps = residua1(n, rvec);	
		
                oldeps = eps;

                printf("eps: %e\n", oldeps);

		if(eps <= EPS){
			count++;
			break;
		}else if(isinf(eps)){
			printf("eps: %e", oldeps);
        		break;
      		}else if (isnan(eps)){
			printf("eps: %e", oldeps);
      			break;
      		}

		//for(i=0, rr=0, pp=0; i<n; i++){	
		//	rr += rvec[i] * rvec[i];
		//	pp += pvec[i] * pvec[i];
		//}
		//beta = rr/pp;

		beta = in(n, rvec,rvec) / in(n, pvec,pvec);

		for(i=0; i<n; i++){
			pvec[i] = rvec[i] + beta * pvec[i];
		}
		count++;
	}
	printf("count: %d\n   eps: %e\n", count, eps);
	free(pvec);free(rvec);free(apvec);
}

int main(int ac, char *av[]){

   int i,j,k;
   int n;
   double *adim;
   double *bvec;
   double *xvec = new double [n]();
   
   char data[30] = "/share/data/matrix/mt10";

   if(ac > 1){
     strcpy(data, av[1]);
   }

   printf("%s\n", data);

   readmats(data, n, adim, bvec);

   for(k=0; k<n; k++){
	for(j=0; j<n; j++){
           //printf("xvec[%2d][%2d] = %lf", k, k, adim[k*n+j], adim[k*n+j]);
   	   printf("%lf	", adim[k*n+j]);
	}
	puts("");
   }

   //cg(n, adim, bvec, xvec);
   gs(n, adim, bvec, xvec);

   for(k=0; k<n; k++){
       printf("xvec[%2d] = %lf\n", k, xvec[k]);
   }

   free(adim);free(bvec);free(xvec);

   return 0;
}
