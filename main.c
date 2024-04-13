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
    int n=1;
    initialiser_parametres(); //initialisations des parametres

    float* x0 = (float*)malloc(Nx*Ny* sizeof(float));
    for (int i=0 ;i<Nx*Ny ;i++)
    {
        x0[i]=1.;
    }
    for (int i=0 ;i<Nx*Ny ;i++)
    {
        printf("%f\n",produit_MV(x0)[i]);
    }

    float* vect_u = (float*)malloc(Nx*Ny* sizeof(float));
    for (int i=0 ;i<Nx*Ny ;i++)
    {
        vect_u[i]=u0(maillage(i,0),maillage(i,1));
    }
    float* vect_un = (float*)malloc(Nx*Ny* sizeof(float));

    switch(time_scheme) {


                    case 1://Explicit
                        double max_vx = abs(v(xmin, ymin)[0]); // Initialisation avec le coin en bas à gauche
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
                        double max_vy = abs(v(xmin, ymin)[1]); // Initialisation avec le coin en bas à gauche
                        for (int i = 0; i < Nx; i++) {
                            double x = xmin + i * dx;
                            for (int j = 0; j < Ny; j++) {
                                double y = ymin + j * dy;
                                double current_vy = abs(v(x, y)[1]);
                                if (current_vy > max_vx) {
                                    max_vy = current_vy;
                                }
                            }
                        }

                        CFL=(max_vx/dx)+(max_vy/dy);
                        float dt=1/CFL;

                        double T=0;
                        while(T<Tf){
                            
                            copierTableau(produit_MV(vect_u),vect_u1,Nx*Ny);
                            scalaireMultiplieTableau(dt,vect_un,vect_un,Nx*Ny);
                            sommeTableaux(vect_u,vect_un,vect_un,Nx*Ny);
                            T+=dt;
                        }

                        break; 


                    case 2://Implicit


                        break;

        
        }


}


