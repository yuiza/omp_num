#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "./header/init.h"

#define NMAX 100
#define EPS 1.0e-10
#define N 3


void gs(int n, double **adim, double *bvec, double *xvec){
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
               xvec[k] -= adim[k][i] * xvec[i];
         }   
         xvec[k] =  xvec[k] / adim[k][k];   

         //εチェック
         eps += fabs(tmp - xvec[k]);
      }
      printf("%d %e   \n", count, eps);
      if(eps < EPS){
         count++;
         // printf("\nEND  %e ", eps);
         break;
      }else if(isinf(eps)){
         // printf("\nERROR : eps infinity\n");
         break;
      }
      count++;
   }
   //printf("%d\n", count);
}

int main(int ac, char *av[]){

   int i,j,k;
   int n = N;
   double **admin;
   double *bvec;
   double *xvec;

   if(ac > 1){
      n = atoi(av[1]);
      // printf("n=%d\n", n);
   }

   admin = (double **)malloc(sizeof(double*)*n);
   for(i=0; i<n; i++){
      admin[i] = (double *)malloc(sizeof(double)*n);
   }
   bvec  = (double *)malloc(sizeof(double)*n);
   xvec  = (double *)malloc(sizeof(double)*n);

   //init A[][], b[], x[]
   init(n, admin, bvec, xvec);

   gs(n, admin, bvec, xvec);

   // for(k=0; k<n; k++){
   //    printf("xvec[%2d] = %e = %lf\n", k,xvec[k], xvec[k]);
   // }
   // printf("xvec[%2d] = %e = %lf\n", n-1,xvec[n-1], xvec[n-1]);


   free(admin);free(bvec);free(xvec);

   return 0;
}