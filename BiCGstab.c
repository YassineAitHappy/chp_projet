#include "BiCGstab.h"
#include "algebre.h"
#include "fonction.h" 
#include "math.h"
#include "stdio.h"
#include "stdlib.h"



void bicgstab(double *b, double *x0, int N, double tol, int max) {
    double *r = (double *)malloc(N * sizeof(double));
    double *r_hat = (double *)malloc(N * sizeof(double));
    double *p = (double *)malloc(N * sizeof(double));
    double *v = (double *)malloc(N * sizeof(double));
    double *s = (double *)malloc(N * sizeof(double));
    double *t = (double *)malloc(N * sizeof(double));
    double *h = (double *)malloc(N * sizeof(double));
    double rho, alpha, omega, beta;

    // Initial values
    copierTableau(produit_MV(x0), r,N); // r = A*x0
    differenceTableaux(b, r, r, N); // r = b - r
    copierTableau(r, r_hat, N); // r_hat = r
    copierTableau(r, p, N); // p = r
    rho = produitScalaire(r_hat, r, N); // rho = (r_hat, r)

    for (int i = 0; i < max; i++) {
        // 1. v = A*p_i-1
        copierTableau(produit_MV(p), v,N);
        
        // 2. alpha = rho_i-1 / (r_hat, v)
        alpha = rho / produitScalaire(r_hat, v, N);
        
        // 3. h = x_i-1 + alpha*p_i-1
        scalaireMultiplieTableau(alpha, p, h, N); // h = alpha*p
        sommeTableaux(x0, h, h, N); // h = x0 + h
        
        // 4. s = r_i-1 - alpha*v
        scalaireMultiplieTableau(-alpha, v, s, N); // s = -alpha*v
        sommeTableaux(r, s, s, N); // s = r + s
        
        // Check for convergence
        if (norm(s, N) < tol) {
            copierTableau(h, x0, N);
            break;
        }

        // 6. t = A*s
        copierTableau(produit_MV(s), t,N);
        
        // 7. omega = (t,s) / (t,t)
        omega = produitScalaire(t, s, N) / produitScalaire(t, t, N);
        
        // 8. x_i = h + omega*s
        scalaireMultiplieTableau(omega, s, x0, N); // x0 = omega*s
        sommeTableaux(h, x0, x0, N); // x0 = h + x0
        
        // 9. r_i = s - omega*t
        scalaireMultiplieTableau(-omega, t, r, N); // r = -omega*t
        sommeTableaux(s, r, r, N); // r = s + r
        
        // Check for convergence
        if (norm(r, N) < tol) {
            break;
        }

        // 10. rho_i = (r_hat, r_i)
        double rho_new = produitScalaire(r_hat, r, N);
        
        // 11. beta = (rho_i/rho_i-1)(alpha/omega)
        beta = (rho_new / rho) * (alpha / omega);
        
        // 12. p_i = r_i + beta(p_i-1 - omega*v)
        scalaireMultiplieTableau(-omega, v, v, N); // v = -omega*v
        sommeTableaux(p, v, p, N); // p = p + v
        scalaireMultiplieTableau(beta, p, p, N); // p = beta*p
        sommeTableaux(r, p, p, N); // p = r + p
        
        rho = rho_new;
    }

    free(r);
    free(r_hat);
    free(p);
    free(v);
    free(s);
    free(t);
    free(h);
}
