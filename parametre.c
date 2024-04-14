#include "parametre.h"
#include <stdio.h>
 
int space_scheme=2, time_scheme=1;
double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
double dt_imp=0.01;
int Nx=100, Ny=100, cas=1;


void initialiser_parametres() {
    switch (cas) {
        case 1:
            xmin=-1.0;
            xmax=1.0;
            ymin=-0.5;
            ymax=0.5;
            Tf=2;
            dx=(xmax-xmin)/Nx;
            dy=(ymax-ymin)/Ny;
            break;
        case 2:
            xmin=-1;
            xmax=1;
            ymin=-1;
            ymax=1;
            Tf=4;
            dx=(xmax-xmin)/Nx;
            dy=(ymax-ymin)/Ny;
            break;
        case 3:
            xmin=-1;
            xmax=1;
            ymin=-1;
            ymax=1;
            Tf=5;
            dx=(xmax-xmin)/Nx;
            dy=(ymax-ymin)/Ny;
            break;
        default:
            printf("choisir un cas!");
            break;
    }
}