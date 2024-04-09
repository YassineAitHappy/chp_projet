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




switch(time_scheme) {


                case 1://Explicit
           
                    break; 


                case 2://Implicit


                    break;

     
     }











   
    }


