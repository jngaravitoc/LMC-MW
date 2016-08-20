#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_legendre.h>


double Anl_tilde(int(n),int(l)){
    double K_nl, factor, A_nl;
    double gamma_factor;

    K_nl = 0.5*n*(n+4.*l+3.) + (l+1.)*(2.*l+1.);
    factor = pow(2,8.*l+6.) / (4.*M_PI*K_nl);
    gamma_factor = pow(gsl_sf_gamma(2.0*l+1.5),2)/gsl_sf_gamma(n+4.*l+3.);
    A_nl =-factor* gsl_sf_fact(n)*(n+2.*l+1.5)*gamma_factor/(4.0*M_PI*K_nl);
    //printf("%f \n", A_nl);
    return A_nl;
}

double phi_nl_f(double(r), int(n), int(l)){
    double factor, s, C_n;
    factor = pow(r, l)*pow((1+r), (-2*l-1))*pow(4*M_PI, 0.5);
    s = (r-1)/(r+1);
    C_n = gsl_sf_gegenpoly_n(n,2.*l+0.5, s);
    //printf("%f \n", C_n);
    return -factor*C_n;
}

double phi_nlm_f(double(r), double(theta),int(n), int(l), int(m)){
    double Y_lm, phi_nl;
    Y_lm = gsl_sf_legendre_sphPlm(l, m, cos(theta));
    phi_nl = phi_nl_f(r, n, l);
    //printf("%f \n", Y_lm);
    return phi_nl*Y_lm;
}

void cov_matrix(double(x), double(y), double(z), double(M), \
               int(nmax), int(lmax), double(r_s)){

    int n, l, m, dm0;
    double r, theta, phi;
    double phi_nlm, A_nl;
    double All_angular;
    r = pow(pow(x, 2) + pow(y, 2) + pow(z,2),0.5);
    theta = acos(z/(r*r_s));
    phi = atan2(y, x);
    //Anl_tilde(2,2);
    phi_nl_f(10, 2, 2);
    for(n=0;n<=nmax;n++){
        for(l=0;l<=lmax;l++){
            for(m=0;m<=l;m++){
            //printf("m = %d, l=%d , cos_theta=%f\n", m, l, cos(theta));
            phi_nlm = phi_nlm_f(r, theta, n, l, m);
            printf("%f \n", phi_nlm);
            if(m==0){
            dm0=0;
            }
            else{
            dm0=1.0;
            }
            A_nl = Anl_tilde(n,l);
            All_angular = phi_nlm*M*cos(m*phi);
            S = (2-dm0)*A_nl*All_angular;
            printf("%f \n", S);
            }
        }
    }
}

int main(){
     cov_matrix(10.,10., 10, 1, 1, 1, 1);
}
