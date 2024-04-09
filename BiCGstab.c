#include "BiCGstab.h"


void bicgstab(float **A, float *b, float *x, int N, float tol, int max_iter) {
    float *r = (float *)malloc(N * sizeof(float));
    float *r_hat = (float *)malloc(N * sizeof(float));
    float *v = (float *)malloc(N * sizeof(float));
    float *p = (float *)malloc(N * sizeof(float));
    float *s = (float *)malloc(N * sizeof(float));
    float *t = (float *)malloc(N * sizeof(float));
    float rho = 1.0, alpha = 1.0, omega = 1.0, rho_1, beta;
    float norm_b = 0.0, norm_r = 0.0, tol_sq = tol * tol;
    int i, iter = 0;

    // Initialisation
    produitMatriceVecteur(A, x, r, N); // r = A*x
    for (i = 0; i < N; ++i) {
        r[i] = b[i] - r[i]; // r = b - A*x
        r_hat[i] = r[i];
        norm_b += b[i] * b[i];
    }
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
            beta = (rho / rho_1) * (alpha / omega);
            for (i = 0; i < N; ++i) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
        }

        produitMatriceVecteur(A, p, v, N); // v = A*p
        alpha = rho / produitScalaire(r_hat, v, N);

        for (i = 0; i < N; ++i) s[i] = r[i] - alpha * v[i]; // s = r - alpha*v

        produitMatriceVecteur(A, s, t, N); // t = A*s
        omega = produitScalaire(t, s, N) / produitScalaire(t, t, N);

        for (i = 0; i < N; ++i) {
            x[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }

        norm_r = produitScalaire(r, r, N);
        if (sqrt(norm_r) < tol) break; // Convergence check
    }

    if (iter >= max_iter) {
        printf("BiCGSTAB did not converge within the maximum number of iterations\n");
    }

    free(r); free(r_hat); free(v); free(p); free(s); free(t);
}
