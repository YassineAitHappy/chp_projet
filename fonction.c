#include "fonction.h"
#include "math.h"
#include "parametre.h"
#include <stdlib.h> 
#include <stdio.h>

extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern int Nx, Ny, cas;


double u0(double x, double y){ 
    double result;
    switch (cas) {
        case 1://translation gaussienne
            result=exp(-(pow(x,2)+pow(y,2))/0.0075);
            break;
        case 2://translation cylindre
            if (sqrt(x*x+y*y) < 0.4){
                result = 1.;
            }
            else{
                result = 0.;
            }
            break;
        case 3://rotation cylindre
            if (sqrt(x*x+y*y) < 0.4){
                result = 1.;
            }
            else{
                result = 0.;
            }
            break;
        default:
            printf("choisir un cas!");
            break;
    }
    return result;
    }

double* v(double x, double y){
    double* result = (double*)malloc(2 * sizeof(double));
    switch (cas) {
        case 1://translation gaussienne
            result[0]=1;
            result[1]=0;
            break;
        case 2://translation cylindre
            result[0]=0.5;
            result[1]=0.5;
            break;
        case 3://rotation cylindre
            result[0] = -y;
            result[1] = x;
            break;
        default:
            printf("choisir un cas!");
            break;
    }
    return result;
}




double alpha( double x,double y){
    double res;
    // res=1.0;
    res = v(x,y)[0]/(2*dx);
    return res;
}
double beta( double x,double y){
    double res;
    //res=2.0;
    res = v(x,y)[1]/(2*dy);
    return res;
}
double gama( double x,double y){
    double res;
    //res=3.0;
    res = -alpha(x,y)-beta(x,y);
    return res;
}

int modulo_pos(int k,int d){
    int r;
    r=k%d;
    if (k<0){
        r=r+d;

    }
    return r;
}