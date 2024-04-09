#include "parametre.h"

 
 int space_scheme, time_scheme;
double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
 int Nx=4, Ny=3, cas=1;


void initialiser_parametres() { 
    switch (cas) {
        case 1:
            xmin=-1;
            xmax=1;
            ymin=-1;
            ymax=1;
            Tf=4;
            break;
        case 2:
            xmin=-1;
            xmax=1;
            ymin=-1;
            ymax=1;
            Tf=5;
            break;
        default:
            xmin=-1;
            xmax=1;
            ymin=-0.5;
            ymax=0.5;
            Tf=2;
            break;
    }
}

