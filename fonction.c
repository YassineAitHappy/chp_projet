#include "fonction.h"
#include "math.h"
#include "parametre.h"
#include <stdlib.h> 

extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern int Nx, Ny, cas;

float u0(float x, float y){ 
    float result;
    switch (cas) {
        case 1://translation cylindre
            if (sqrt(x*x+y*y) < 0.4){
                result = 1.;
            }
            else{
                result = 0.;
            }
            break;
        case 2://rotation cylindre
            if (sqrt(x*x+y*y) < 0.4){
                result = 1.;
            }
            else{
                result = 0.;
            }
            break;
        default://translation gaussienne
            result=exp((-pow(x,2)/0.0075)-pow(y,2)/0.0075);
            break;
    }
    return result;
    }

float* v(float x, float y){
    float* result = (float*)malloc(2 * sizeof(float));
    switch (cas) {
        case 1://translation cylindre
            result[0]=0.5;
            result[1]=0.5;
            break;
        case 2://rotation cylindre
            result[0] = -y;
            result[1] = x;
            break;
        default://translation gaussienne
            result[0]=1;
            result[1]=0;
            break;
    }
    return result;
}



float alpha( float x,float y){
return 1.0;
}
float beta( float x,float y){
    return 2.0;

}
float gama( float x,float y){
    return 3.0;
}