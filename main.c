#include "stdio.h"
#include "stdlib.h"
#include "fonction.h"
#include "algebre.h"
#include "BiCGstab.h"
#include "parametre.h"
#include "math.h"
extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern int Nx, Ny, cas;
int main(int argc, char* argv[]) {
    initialiser_parametres();

    // testing

    double dt_test=0.01;
    double* x0 = (double*)malloc(Nx*Ny* sizeof(double));//ce vecteur permet l'affichage d'une colonne de la matrice
    for (int i=0 ;i<Nx*Ny ;i++)
    {   
        if (i==3){
            x0[i]=1;
        }
        else{
        x0[i]=0.;
        }
        
    }
    for (int i=0 ;i<Nx*Ny ;i++)
    {   
        if (space_scheme==1){//A=I-dt*B
            //printf("produit implicit: %f\n",x0[i]-dt_test*produit_MV(x0)[i]);
            
        }
        else if(space_scheme==2){//A=I+dt*B
            //printf("produit explicit: %f\n",x0[i]+dt_test*produit_MV(x0)[i]);
           
        }
    }

    double* vect_u = (double*)malloc(Nx*Ny* sizeof(double));
    for (int i=0 ;i<Nx*Ny ;i++)
    {   
        vect_u[i]=u0(maillage(i,0),maillage(i,1));//vecteur initiale
    }
    double* vect_un = (double*)malloc(Nx*Ny* sizeof(double));


    switch(time_scheme) {
    case 1: //Explicit
        {
            //Trouver la cfl
            //trouver V_x max
            double max_vx; // Déclaration de max_vx au début du bloc de code
            max_vx = fabs(v(xmin, ymin)[0]); // Initialisation avec le coin en bas à gauche
            for (int i = 0; i < Nx; i++) {
                double x = xmin + i * dx;
                for (int j = 0; j < Ny; j++) {
                    double y = ymin + j * dy;
                    double current_vx = fabs(v(x, y)[0]);
                    if (current_vx > max_vx) {
                        max_vx = current_vx;
                    }
                }
            }
            //trouver V_y max
            double max_vy; // Déclaration de max_vy au début du bloc de code
            max_vy = fabs(v(xmin, ymin)[1]); // Initialisation avec le coin en bas à gauche
            for (int i = 0; i < Nx; i++) {
                double x = xmin + i * dx;
                for (int j = 0; j < Ny; j++) {
                    double y = ymin + j * dy;
                    double current_vy = fabs(v(x, y)[1]);
                    if (current_vy > max_vy) {
                        max_vy = current_vy;
                    }
                }
            }

            CFL = (max_vx / dx) + (max_vy / dy);
            printf("cfl=%f",CFL);
            double dt = 1.0 / CFL;
            printf("dt=%f\n",dt);
            

            char filenam[120]; // Définissez la taille en fonction du format du nom du fichier
            sprintf(filenam, "./sol.%d.dat", 0); // Formatage du nom du fichier
            FILE* file_explicit = fopen(filenam, "w"); // Ouvrir le fichier en écriture
            if (file_explicit == NULL) {
                printf("Erreur lors de la création du fichier %s\n", filenam);
            return 1;
            }
            //stocker le vecteur initial
            for (int i = 0; i < Nx * Ny; i++) {
                fprintf(file_explicit, "%f %f %f\n", maillage(i, 0), maillage(i, 1), vect_u[i]);
            }
            fprintf(file_explicit, "\n"); // Saut de ligne pour passer à l'itération de temps suivante
            fclose(file_explicit); // Fermer le fichier


            // Boucle temporelle pour l'explicite
            double T_explicit = 0;
            int Nt=0;
            while (T_explicit < Tf) {
                // Mise à jour du vecteur de résolution pour la prochaine étape de temps
                copierTableau(produit_MV(vect_u), vect_un, Nx * Ny);
                scalaireMultiplieTableau(dt, vect_un, vect_un, Nx * Ny);
                sommeTableaux(vect_u, vect_un, vect_u, Nx * Ny);
                T_explicit += dt; // Incrémentation du temps
                Nt+=1;
                char filename[120]; // Définissez la taille en fonction du format du nom du fichier
                sprintf(filename, "./sol.%d.dat", Nt); // Formatage du nom du fichier
                FILE* file_explicit = fopen(filename, "w"); // Ouvrir le fichier en écriture
                if (file_explicit == NULL) {
                    printf("Erreur lors de la création du fichier %s\n", filename);
                return 1;
                }
                for (int i = 0; i < Nx * Ny; i++) {
                    fprintf(file_explicit, "%f %f %f\n", maillage(i, 0), maillage(i, 1), vect_u[i]);
                }
                fprintf(file_explicit, "\n"); // Saut de ligne pour passer à l'itération de temps suivante
                fclose(file_explicit); // Fermer le fichier
            }

            
        }
        break;
        
    case 2: //Implicit
    {
        char filenam[120]; // Définissez la taille en fonction du format du nom du fichier
            sprintf(filenam, "./sol.%d.dat", 0); // Formatage du nom du fichier
            FILE* file_implicit = fopen(filenam, "w"); // Ouvrir le fichier en écriture
            if (file_implicit == NULL) {
                printf("Erreur lors de la création du fichier %s\n", filenam);
            return 1;
    }
            for (int i = 0; i < Nx * Ny; i++) {
                    fprintf(file_implicit, "%f %f %f\n", maillage(i, 0), maillage(i, 1), vect_u[i]);
                }
                fprintf(file_implicit, "\n"); // Saut de ligne pour passer à l'itération de temps suivante
            fclose(file_implicit); // Fermer le fichier
            // Boucle temporelle pour l'explicite
            double T_implicit = 0;
            int Nt=0;
            while (T_implicit < Tf) {
                // Résolution du système d'équations pour la prochaine étape de temps
                bicgstab(vect_u, vect_un, Nx * Ny, 1E-6, 5000);
                copierTableau(vect_un,vect_u,Nx*Ny);

                T_implicit += dt_imp; // Incrémentation du temps
                Nt+=1;
                char filename[120]; // Définissez la taille en fonction du format du nom du fichier
                sprintf(filename, "./sol.%d.dat", Nt); // Formatage du nom du fichier
                FILE* file_implicit = fopen(filename, "w"); // Ouvrir le fichier en écriture
                if (file_implicit == NULL) {
                    printf("Erreur lors de la création du fichier %s\n", filename);
                return 1;
                }
                for (int i = 0; i < Nx * Ny; i++) {
                    fprintf(file_implicit, "%f %f %f\n", maillage(i, 0), maillage(i, 1), vect_u[i]);
                }
                fprintf(file_implicit, "\n"); // Saut de ligne pour passer à l'itération de temps suivante
                fclose(file_implicit); // Fermer le fichier
            }   
            }
        break;
    }

    // Libération de la mémoire allouée
    free(x0);
    free(vect_u);
    free(vect_un);

    return 0;
}


