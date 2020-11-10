#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
#include <iostream>

int main()
{
    int n = 10;
    double *d1 = (double *) calloc(n, sizeof(double));
    double *d2 = (double *) calloc(n-1, sizeof(double));
    double *M = (double *) calloc(n*n, sizeof(double));
    double *M2 = (double *) calloc(n*n, sizeof(double));


    double *d13_r_eigva = (double *) calloc(n, sizeof(double));
    double *d13_i_eigva = (double *) calloc(n, sizeof(double));
    double *d13_eigvec = (double *) calloc(n*n, sizeof(double));

    int info1, info2, info3;

    FILE *fp;

    //for the DSTEV sub-routine-----------------------------------------------

    for(int i = 0; i<n-1; i++)
    {
        d1[i] = 2;
        d2[i] = -1;
    }
    d1[n-1] = 2; //that's out of the loop to avoid a if

    info1 = LAPACKE_dstev(LAPACK_ROW_MAJOR, 'V', n, d1, d2, M, n);

    fp = fopen("VEC_DSTEV.dat", "w");
        for(int j = 0;j<n; j++)
        {
            fprintf(fp,"%lf ", d1[j]); //first line are the eigenvalues
        }
        
        fprintf(fp,"\n \n \n");

        for (int i=0; i<n; i++) 
        {
            for(int j = 0;j<n; j++)
            {
                fprintf(fp,"%lf ", M[i*n+j]); //eigenvectors sorted as columns
            }
            fprintf(fp,"\n");
        }
    fclose(fp);

    free(d1);
    free(d2);
    free(M);
    //-----------------------DSTEV sub-routine-----------------------------------------------

    //for the DSYEV sub-routine-----------------------------------------------------------------------
    for(int i = 0; i<n; i++)
    {
        M[i+n*i] = 2;
    }

    for(int i =0; i<(n-1); i++) //to avoi a if 
    {
        M[i+n*(i+1)] = -1*(n!=(n-1)); //I'll fill only the bottom part
    }

    info2 = LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','L', n, M, n, d1); 

    for(int j = 0;j<n; j++)
    {
        fprintf(fp,"%lf ", d1[j]); //first line are the eigenvalues
    }
    
    fprintf(fp,"\n \n \n");

    for (int i = 0; i<n; i++) 
    {
        for(int j = 0;j<n; j++)
        {
            fprintf(fp,"%lf ",M[i*n+j]); //eigenvectors sorted in columns
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    free(fp);
    free(M);
    free(d2);
    //----------------DSYEV sub-routine-----------------------------------------------------------------------

    //for the DGEEV sub-routine-----------------------------------------------------------------------
    for(int i = 0; i<n; i++)
    {
        M2[i+n*i] = 2;
    }
    for(int i = 0; i<(n-1); i++)
    {
        M[i+n*(i+1)] = -1; //I'll fill only the bottom part
        M[i+1+n*i] = -1;
    }

    info3 = LAPACKE_dgeev(LAPACK_ROW_MAJOR,'N','V', n, M2, n, d13_r_eigva, d13_i_eigva,NULL, n, d13_eigvec, n);

    fp = fopen("VEC_DGEEV.dat", "w");
    for(int j = 0;j<n; j++)
    {
        fprintf(fp,"%lf ", d13_r_eigva[j]); //first line are the eigenvalues
    }

    fprintf(fp,"\n \n \n");
    
    for (int i = 0; i<n; i++) 
    {
        for(int j = 0;j<n; j++)
        {
            fprintf(fp,"%lf ",d13_eigvec[i*n+j]); //eigenvectors sorted in columns
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    free(fp);
    free(M);
    free(M2);
    free(d13_eigvec);
    free(d13_r_eigva);
    free(d13_i_eigva);
    //----------------DGEEV sub-routine-----------------------------------------------------------------------


    return 0;

}