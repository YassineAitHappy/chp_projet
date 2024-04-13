#include "stdio.h"
#include "stdlib.h"
#include "fonction.h"
#include "algebre.h"
#include "BiCGstab.h"
#include "parametre.h"
extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern int Nx, Ny, cas;
int main(int argc, char* argv[]) {
    initialiser_parametres();
    double dt_test=0.01;
    double* x0 = (double*)malloc(Nx*Ny* sizeof(double));
    for (int i=0 ;i<Nx*Ny ;i++)
    {
        x0[i]=1.;
    }
    for (int i=0 ;i<Nx*Ny ;i++)
    {
        printf("%f\n",x0[i]+dt_test*produit_MV(x0)[i]);
    }

    double* vect_u = (double*)malloc(Nx*Ny* sizeof(double));
    for (int i=0 ;i<Nx*Ny ;i++)
    {
        vect_u[i]=u0(maillage(i,0),maillage(i,1));
        printf("xi=%f,uj=%f: %f\n",maillage(i,0),maillage(i,1),vect_u[i]);
    }
    double* vect_un = (double*)malloc(Nx*Ny* sizeof(double));

    switch(time_scheme) {
    case 1: //Explicit
        {
            double max_vx; // Déclaration de max_vx au début du bloc de code
            max_vx = abs(v(xmin, ymin)[0]); // Initialisation avec le coin en bas à gauche
            for (int i = 0; i < Nx; i++) {
                double x = xmin + i * dx;
                for (int j = 0; j < Ny; j++) {
                    double y = ymin + j * dy;
                    double current_vx = abs(v(x, y)[0]);
                    if (current_vx > max_vx) {
                        max_vx = current_vx;
                    }
                }
            }

            double max_vy; // Déclaration de max_vy au début du bloc de code
            max_vy = abs(v(xmin, ymin)[1]); // Initialisation avec le coin en bas à gauche
            for (int i = 0; i < Nx; i++) {
                double x = xmin + i * dx;
                for (int j = 0; j < Ny; j++) {
                    double y = ymin + j * dy;
                    double current_vy = abs(v(x, y)[1]);
                    if (current_vy > max_vy) {
                        max_vy = current_vy;
                    }
                }
            }

            CFL = (max_vx / dx) + (max_vy / dy);
            double dt = 1.0 / 10*CFL;
            printf("dt=%f\n",dt);
            

            FILE* file_explicit = fopen("result_explicit.txt", "w"); // Ouvrir le fichier en écriture
                if (file_explicit == NULL) {
                    printf("Erreur lors de la création du fichier result_explicit.txt\n");
                    return 1;
                }

            // Boucle temporelle pour l'explicite
            double T_explicit = 0;
            while (T_explicit < Tf) {
                // Mise à jour du vecteur de résolution pour la prochaine étape de temps
                copierTableau(produit_MV(vect_u), vect_un, Nx * Ny);
                scalaireMultiplieTableau(dt, vect_un, vect_un, Nx * Ny);
                sommeTableaux(vect_u, vect_un, vect_un, Nx * Ny);

                // Écrire les valeurs du vecteur de résolution dans le fichier pour cette itération de temps
                for (int i = 0; i < Nx * Ny; i++) {
                    fprintf(file_explicit, "%f ", vect_u[i]);
                }
                fprintf(file_explicit, "\n"); // Saut de ligne pour passer à l'itération de temps suivante

                T_explicit += dt; // Incrémentation du temps
            }

            fclose(file_explicit); // Fermer le fichier
        }
        break;
        
    case 2: //Implicit
    {
        FILE* file_implicit = fopen("result_implicit.txt", "w"); // Ouvrir le fichier en écriture
                if (file_implicit == NULL) {
                    printf("Erreur lors de la création du fichier result_implicit.txt\n");
                    return 1;
                }

                // Boucle temporelle pour l'implicite
                double T_implicit = 0;
                while (T_implicit < Tf) {
                    // Résolution du système d'équations pour la prochaine étape de temps
                    bicgstab(vect_u, vect_un, Nx * Ny, 0.0001, 20);

                    // Écrire les valeurs du vecteur de résolution dans le fichier pour cette itération de temps
                    for (int i = 0; i < Nx * Ny; i++) {
                        fprintf(file_implicit, "%f ", vect_u[i]);
                    }
                    fprintf(file_implicit, "\n"); // Saut de ligne pour passer à l'itération de temps suivante

                    T_implicit += dt_imp; // Incrémentation du temps
                }

                fclose(file_implicit); // Fermer le fichier
            }
        break;
    }

    // Libération de la mémoire allouée
    free(x0);
    free(vect_u);
    free(vect_un);

    return 0;
}

