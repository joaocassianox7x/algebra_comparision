#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
#include <armadillo>
#include <iostream>
#include <time.h>

//Here i generata the time spent of each function//

using namespace arma;

int main()
{
    double *d1;
    double *d2;
    double *M;
    double *M2;
    double *M3;

    double *d13_r_eigva;
    double *d13_i_eigva;
    double *d13_eigvec;

    int info;

    clock_t t0, t1;

    int n;
    int m =25;

    double *vec_t1 = (double *) calloc(m, sizeof(double));
    double *vec_t2 = (double *) calloc(m, sizeof(double));
    double *vec_t3 = (double *) calloc(m, sizeof(double));
    double *vec_t4 = (double *) calloc(m, sizeof(double));

    for(int n_int = 0; n_int<m; n_int++)
    {
        n = 100+n_int*100;
        printf("N = %d\n",n);
        d1 = (double *) calloc(n, sizeof(double));
        d2 = (double *) calloc(n-1, sizeof(double));
        M = (double *) calloc(n*n, sizeof(double));
        M2 = (double *) calloc(n*n, sizeof(double));
        M3 = (double *) calloc(n*n, sizeof(double));

        d13_r_eigva = (double *) calloc(n, sizeof(double));
        d13_i_eigva = (double *) calloc(n, sizeof(double));
        d13_eigvec = (double *) calloc(n*n, sizeof(double));

        for(int i = 0; i<n-1; i++)
        {
            d1[i] = 2.;
            d2[i] = -1.;
        }
        d1[n-1] = 2.; //that's out of the loop to avoid a if

        t0 = clock();
        info = LAPACKE_dstevd(LAPACK_ROW_MAJOR, 'V', n, d1, d2, M, n);
        t1=clock();
        vec_t1[n_int] = (double)(t1-t0)/(CLOCKS_PER_SEC);

        for(int i = 0; i<n; i++)
        {
            M2[i+n*i] = 2;
        }

        for(int i =0; i<(n-1); i++) //to avoi a if 
        {
            M2[i+n*(i+1)] = -1*(n!=(n-1)); //I'll fill only the bottom part
        }

        t0 = clock();
        info = LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','L', n, M2, n, d1); 
        t1=clock();
        vec_t2[n_int] = (double)(t1-t0)/(CLOCKS_PER_SEC);
        
        for(int i = 0; i<n; i++)
        {
            M3[i+n*i] = 2;
        }
        for(int i = 0; i<(n-1); i++)
        {
            M3[i+n*(i+1)] = -1;
            M3[i+1+n*i] = -1;
        }
        

        t0 = clock();
        info = LAPACKE_dgeev(LAPACK_ROW_MAJOR,'N','V', n, M3, n, d13_r_eigva, d13_i_eigva,NULL, n, d13_eigvec, n);
        t1=clock();
        vec_t3[n_int] = (double)(t1-t0)/(CLOCKS_PER_SEC);

        
        mat MA(n,n);
        MA.diag(-1).fill(-1.);
        MA.diag(+1).fill( -1.);
        MA.diag(0).fill(2.);

        vec eigval;
        mat eigvec;

        t0 = clock();
        eig_sym(eigval,eigvec, MA);
        t1=clock();
        vec_t4[n_int] = (double)(t1-t0)/(CLOCKS_PER_SEC);

        free(d1);
        free(d2);
        free(M);
        free(M2);
        free(M3);
        free(d13_r_eigva);
        free(d13_i_eigva);
        free(d13_eigvec);

    }

    FILE *fp;

    fp = fopen("tempos.dat","w");
    for(int i=0; i<m;i++)
    {
        fprintf(fp,"%lf %lf %lf %lf\n", vec_t1[i], vec_t2[i], vec_t3[i],vec_t4[i]);
    }
    fclose(fp);

    printf("%d",info);
    return 0;

}