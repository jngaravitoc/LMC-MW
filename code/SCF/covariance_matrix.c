/*
Author: J. Nicolas Garavito-Camargo
08/19/2016
University of Arizona.

Code to compute the covariance matrix of the
coefficients of the basis expansions functions.

usage:

the code reads an input file

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_legendre.h>

double Anl_tilde(int n, int l);
double phi_nl_f(double r, int n, int l);
double phi_nlm_f(double r, double theta, int n, int l, int m);
double sum_angular(int n_points, double *r, double *theta, double *phi, double *M, int n, int l, int m);
void read_data(char *filename, int n_points, double *r, double *theta, double *phi, double *m, double r_s);
void cov_matrix(int n_points, double *r , double *theta , double *phi, double *m, int max, int lmax);


void coefficients(int n_points, double *r , double *theta , double *phi, double *m, int max, int lmax);

int main(void){
     double *r=NULL;
     double *theta=NULL;
     double *phi=NULL;
     double *M=NULL;

     /* Global variables */
     /* To do: put this as input parameters */

     int n_points=425584;
     int nmax=10;
     int lmax=1;

     char filename[20]="spherical_halo.txt";
     double r_s = 40.85;

     // ------------------------

     printf("%s \n", filename);

     /* Allocating memory for pointers */
     r = malloc(n_points*sizeof(long double));
     theta = malloc(n_points*sizeof(long double));
     phi = malloc(n_points*sizeof(long double));
     M = malloc(n_points*sizeof(long double));

     read_data(filename, n_points, r, theta, phi, M, r_s);
     cov_matrix(n_points, r, theta, phi, M, nmax, lmax);;
     return 0;
}

double Anl_tilde(int n ,int l){
    double K_nl, factor, A_nl;
    double gamma_factor;

    K_nl = 0.5*n*(n+4.*l+3.) + (l+1.)*(2.*l+1.);
    factor = pow(2,8.*l+6.) / (4.*M_PI*K_nl);
    gamma_factor = pow(gsl_sf_gamma(2.0*l+1.5),2)/gsl_sf_gamma(n+4.*l+3.);
    A_nl =-factor* gsl_sf_fact(n)*(n+2.*l+1.5)*gamma_factor;
    return A_nl;
}

double phi_nl_f(double r, int n, int l){
    double factor, s, C_n;
    factor = pow(r, l) * pow((1.+r), (-2.*l-1.)) * pow(4*M_PI, 0.5);
    s = (r-1)/(r+1);
    C_n = gsl_sf_gegenpoly_n(n,2*l+1.5, s);
    return -factor*C_n;
}

/* Function that computes the potential see Eq. in Lowing*/
double phi_nlm_f(double r, double theta ,int n, int l, int m){
    double Y_lm, phi_nl;
    Y_lm = gsl_sf_legendre_sphPlm(l, m, cos(theta));
    phi_nl = phi_nl_f(r, n, l);
    return phi_nl*Y_lm;
}

/* function that sums the angular terms over all the particles */
double sum_angular(int n_points, double *r, double *theta, double *phi, double *M, int n, int l, int m){
    double all_angular=0;
    double phi_nlm;
    int i;
    for(i=0;i<=n_points;i++){
    phi_nlm = phi_nlm_f(r[i], theta[i], n, l, m);
    //printf("%f \n", phi_nlm);
    all_angular += phi_nlm*cos(m*phi[i])*M[i];
    }
    return all_angular;
}

void cov_matrix(int n_points, double *r , double *theta , double *phi, double *M, int nmax, int lmax){

    int n, l, m, dm0;
    double A_nl;
    double All_angular;
    double S;
    for(n=0;n<=nmax;n++){
        for(l=0;l<=lmax;l++){
            for(m=0;m<=l;m++){

            if(m==0){
            dm0=1.0;
            }
            else{
            dm0=0.0;
            }

            A_nl = Anl_tilde(n,l);
            All_angular = sum_angular(n_points, r, theta, phi, M, n, l, m);
            S = (2-dm0)*A_nl*All_angular;

            printf("%f \t \n", S);
            }
        }
    }
}


void coefficients(int n_points, double *r , double *theta , double *phi, double *M, int nmax, int lmax){

    int n, l, m, dm0;
    double A_nl;
    double All_angular;
    double S;
    for(n=0;n<=nmax;n++){
        for(l=0;l<=lmax;l++){
            for(m=0;m<=l;m++){

            if(m==0){
            dm0=1.0;
            }
            else{
            dm0=0.0;
            }

            A_nl = Anl_tilde(n,l);
            All_angular = sum_angular(n_points, r, theta, phi, M, n, l, m);
            S = (2-dm0)*A_nl*All_angular;

            //intf("%f \t \n", S);
            }
        }
    }
}
void read_data(char *filename, int n_points, double *r, double *theta, \
               double *phi, double *M, double r_s){

    FILE *in;
    double X, Y, Z, m;
    int i;

    in = fopen(filename, "r");

    /* Checking if the file is opening*/
    if(!in){
        printf("Problem opening file %s \n", filename);
        exit(1);
    }

    for(i=0;i<=n_points;i++){
        fscanf(in, "%lf %lf %lf %lf \n", &X, &Y, &Z, &m);
        r[i] = pow(pow(X, 2) + pow(Y, 2) + pow(Z,2),0.5)/r_s;
        theta[i] = acos(Z/(r[i]*r_s));
        phi[i] = atan2(Y, X);
        M[i] = m;
    }
    fclose(in);
}
