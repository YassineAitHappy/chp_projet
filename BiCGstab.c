#include "BiCGstab.h"
#include "algebre.h"
#include "fonction.h" 
#include "math.h"
#include "stdio.h"
#include "stdlib.h"


void bicgstab( double *b, double *x, int N, double tol, int max_iter) {
    double *r = (double *)malloc(N * sizeof(double));
    double *r_hat = (double *)malloc(N * sizeof(double));
    double *v = (double *)malloc(N * sizeof(double));
    double *p = (double *)malloc(N * sizeof(double));
    double *s = (double *)malloc(N * sizeof(double));
    double *t = (double *)malloc(N * sizeof(double));
    double rho = 1.0, coeff1 = 1.0, omega = 1.0, rho_1,coeff2 ;
    double norm_b = 0.0, norm_r = 0.0, tol_sq = tol * tol;
    int i, iter = 0;

    // Initialisation
    differenceTableaux(b,produit_MV(x),r,N); // r = b-A*x
    copierTableau(r,r_hat,N);
    norm_b=norm( b, N);
    if (sqrt(norm_b) < tol) {
        free(r); free(r_hat); free(v); free(p); free(s); free(t);
        return;
    }

    for (iter = 0; iter < max_iter; ++iter) {
        rho_1 = rho;
        rho = produitScalaire(r_hat, r, N);

        if (fabs(rho) < tol) {
            printf("BiCGSTAB failed to converge\n");
            break;
        }

        if (iter == 0) {
            copierTableau(r, p, N);
        } else {
            coeff2 = (rho / rho_1) * (coeff1 / omega);
            for (i = 0; i < N; ++i) {
                p[i] = r[i] + coeff2 * (p[i] - omega * v[i]);
            }
        }
        
        copierTableau(produit_MV(p),v,N); // v = A*p
        coeff1 = rho / produitScalaire(r_hat, v, N);

        scalaireMultiplieTableau(coeff1,v,s,N);
        differenceTableaux(r,s,s,N);// s = r - coeff1*v

        copierTableau(produit_MV(s),t,N); // t = A*s
        omega = produitScalaire(t, s, N) / produitScalaire(t, t, N);

        scalaireMultiplieTableau(coeff1,p,p,N);
        scalaireMultiplieTableau(omega,s,s,N);
        sommeTableaux(p,s,x,N);//x=coeff1*p+omega*s
        
        scalaireMultiplieTableau(omega,t,r,N);
        differenceTableaux(s,r,r,N);;//r=s-omega*t

        norm_r =norm(r,N);
        if (sqrt(norm_r) < tol) break; // Convergence check
    }

    if (iter >= max_iter) {
        printf("BiCGSTAB did not converge within the maximum number of iterations\n");
    }

    free(r); free(r_hat); free(v); free(p); free(s); free(t);
}