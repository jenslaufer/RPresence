#include <stdio.h>
#include <math.h>
#include "SingAnal.h"
bool InvertMatrix2(double **a, int n)
/* Inverts matrix using Cholesky decomposition*/
{ FILE *g;
    int ii,jj,kk,i,j;    double sum, *p;    p = new double[n]; bool rv;
    extern TSingSpec SingSpec;    g=SingSpec.g;
    for (ii=0; ii<n; ii++) {
        for (jj=ii; jj<n; jj++) {
            for (sum=a[ii][jj],kk=ii-1;kk>=0;kk--) { sum -= a[ii][kk]*a[jj][kk];  }
            if (ii == jj) {
                if (sum==0.0 || (sum<0.0 && sum>-.1e-7)) sum=.1e-20;
                if (sum <= 0.0) {
                    fprintf(g,"\n\n**********************   WARNING   ***************************\n");
                    fprintf(g,"Variance-covariance matrix has not been computed successfully.\n");
                    fprintf(g,"Ignore matrix values in the below output.\n");
					fprintf(g,"    (sum=%f)\n",sum);
                    fprintf(g,"--------------------------------------------------------------\n");
                    delete[] p;
                    return(false);
                }
                p[ii]=sqrt(sum);
            } else {
                a[jj][ii]=sum/p[ii];
            }
        }
    }
    for (ii=0; ii<n; ii++) {
        a[ii][ii]=1.0/p[ii];
        for (jj=ii+1; jj<n; jj++) {
            for (kk=ii,sum=0.; kk<jj; kk++) { sum -= a[jj][kk]*a[kk][ii]; }
            a[jj][ii]=sum/p[jj];  a[ii][jj]=0.0;
        }
    }
    
    double **temp;   temp = new double*[n];
    for (ii=0; ii<n; ii++) {
        temp[ii] = new double[n];
        for (jj=0; jj<n; jj++) {
            temp[ii][jj] = 0.0;
            for (kk=0; kk<n; kk++) {temp[ii][jj] += a[kk][ii]*a[kk][jj]; }
        }
    }
    for (ii=0; ii<n; ii++) {
        for (jj=0; jj<n; jj++) { a[ii][jj] = temp[ii][jj]; }
        delete[] temp[ii];
    }
    if(SingSpec.Verbose>1) fprintf(g,"mat*inverse:\n"); sum=0;
    
    for (i=0; i<n; i++) { 
        for (j=0; j<n; j++) { 
            if(SingSpec.Verbose>1) fprintf(g," %f",a[i][j]); 
            sum+=fabs(a[i][j]);
        }
        if(SingSpec.Verbose>1) fprintf(g,"\n");
    }	   
    rv=(fabs(sum-n)>.001); 
    if (rv) {
        fprintf(g,"*** Matrix Inversion Error (sum=%f) ***\n",sum);	
        fprintf(g,"\n\n**********************   WARNING   ***************************\n");
        fprintf(g,"Variance-covariance matrix has not been computed successfully.\n");
        fprintf(g,"Ignore variances and covariances in the output below.\n");
        fprintf(g,"--------------------------------------------------------------\n");
    }
    delete[] temp;  delete[] p;
    return(rv);
}
int InvertMatrix(double **a, int n)  {  
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
    int i,j,k; double x,sum=0; double **b,**c; bool rv; 
    if (n == 1) return(0);  // must be of dimension >= 2
	b=new double*[n]; c=new double*[n];
    if(SingSpec.Verbose>1) fprintf(g,"inverting the following matrix...\n");
    for (i=0; i<n; i++) { 
        b[i]=new double[n]; c[i]=new double[n];
        for (j=0; j<n; j++) { 
            b[i][j]=a[i][j]; c[i][j]=0;
            if(SingSpec.Verbose>1) fprintf(g," %10.6f",a[i][j]);
        }
        if(SingSpec.Verbose>1) fprintf(g,"\n");
    }
    for (i=1; i < n; i++) a[0][i] /= a[0][0]; // normalize row 0
    for (i=1; i < n; i++)  { 
        for (j=i; j < n; j++)  { // do a column of L
            sum = 0.0;
            for (k = 0; k < i; k++) sum += a[j][k] * a[k][i];
            a[j][i] -= sum;
        }
        if (i == (n-1)) continue;
        for (j=i+1; j < n; j++)  {  // do a row of U
            sum = 0.0;
            for (k = 0; k < i; k++) sum += a[i][k]*a[k][j];
            a[i][j] = (a[i][j]-sum) / a[i][i];
        }
    }
    for ( i = 0; i < n; i++ )  // invert L
        for ( j = i; j < n; j++ )  {
            x = 1.0;
            if ( i != j ) {
                x = 0.0;
                for ( k = i; k < j; k++ ) x -= a[j][k]*a[k][i];
            }
            a[j][i] = x / a[j][j];
        }
    for ( i = 0; i < n; i++ )   // invert U
        for ( j = i; j < n; j++ )  {
            if ( i == j ) continue;
            sum = 0.0;
            for ( k = i; k < j; k++ ) sum += a[k][j]*( (i==k) ? 1.0 : a[i][k] );
            a[i][j] = -sum;
        }
    for ( i = 0; i < n; i++ )   // final inversion
        for ( j = 0; j < n; j++ )  {
            sum = 0.0;
            for ( k = ((i>j)?i:j); k < n; k++ ) sum += ((j==k)?1.0:a[j][k])*a[k][i];
            a[j][i] = sum;
        }
    if(SingSpec.Verbose>1) fprintf(g,"inverse:\n");
    for (i=0; i<n; i++) { 
        for (j=0; j<n; j++) { 
			if(SingSpec.Verbose>1) fprintf(g," %f",a[i][j]);
            sum=0; for (k=0; k<n; k++) sum+=b[i][k]*a[k][j];
            c[i][j]=sum;
        }
        if(SingSpec.Verbose>1) fprintf(g,"\n");
    }
    if(SingSpec.Verbose>1) fprintf(g,"mat*inverse:\n"); sum=x=0;
    for (i=0; i<n; i++) { 
        for (j=0; j<n; j++) { 
            if(SingSpec.Verbose>1) fprintf(g," %f",c[i][j]); 
            sum+=fabs(c[i][j]);
        }
        if(SingSpec.Verbose>1) fprintf(g,"\n");
		delete [] b[i]; delete [] c[i];
    }	
    delete [] b; delete [] c; 
    rv=(fabs(sum-n)>.001);
    if (rv) {
        fprintf(g,"*** Matrix Inversion Error (sum=%f) ***\n",sum);	
        fprintf(g,"\n\n**********************   WARNING   ***************************\n");
        fprintf(g,"Variance-covariance matrix has not been computed successfully.\n");
        fprintf(g,"Ignore variances and covariances in the output below.\n");
        fprintf(g,"--------------------------------------------------------------\n");
    }
    return(rv);
}
