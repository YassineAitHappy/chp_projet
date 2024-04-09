#include "algebre.h"
#include "math.h"
#include "parametre.h"
#include <stdlib.h> 
#include "fonction.h"
extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern int Nx, Ny, cas;

float* produit_MV( float* vecteur) {

    float* resultat = (float*)malloc(Nx*Ny* sizeof(float)); 


    // Calcul du produit matrice-vecteur
    if (vecteur== NULL) {
        
        return NULL;
    }
    else{
     switch(time_scheme) {
        case 1://Explicit
    
            switch(space_scheme) {

                case 2://Upwind

                    for (int I=1;I<Nx*Ny;I++){
                        
                            
                        if (I==0)
                        {
                            resultat[I]=gama(maillage(I,0),maillage(I,1))*vecteur[I]+alpha(maillage(I+Nx-1,0),maillage(I+Nx-1,1))*vecteur[I+Nx-1]+beta(maillage(I+(Ny-1)*Nx,0),maillage(I+(Ny-1)*Nx,1))*vecteur[I+(Ny-1)*Nx];

                        }
                        else if((0<I)&&(I<Nx) )
                        {
                           resultat[I]=gama(maillage(I,0),maillage(I,1))*vecteur[I]+alpha(maillage(I-1,0),maillage(I-1,1))*vecteur[I-1]+beta(maillage(I+(Ny-1)*Nx,0),maillage(I+(Ny-1)*Nx,1))*vecteur[I+(Ny-1)*Nx];
                        }
                        else if(I%Nx==0)
                        {
                            resultat[I]=gama(maillage(I,0),maillage(I,1))*vecteur[I]+alpha(maillage(I-1,0),maillage(I-1,1))*vecteur[I-1]+beta(maillage(I-Nx,0),maillage(I-Nx,1))*vecteur[I-Nx]+alpha(maillage(I+Nx-1,0),maillage(I+Nx-1,1))*vecteur[I+Nx-1];
                        }
                        else
                        {
                            resultat[I]=gama(maillage(I,0),maillage(I,1))*vecteur[I]+alpha(maillage(I-1,0),maillage(I-1,1))*vecteur[I-1]+beta(maillage(I-Nx,0),maillage(I-Nx,1))*vecteur[I-Nx];
                        }

                        }





                                









                    break;
            }
            break;

        case 2://Implicit
            switch(space_scheme) {

                case 1://Centré
           
                    break; 
        
             }
          
            break;
     }

    
   
    }
    return resultat; // Retourne le pointeur vers le vecteur résultat
}


int indexe(int i,int j){ // focntion qui donne l'indice de l'element stocké en fonction de i et j
    return (j-1)*Nx+i-1;
} 
int couple(int I,int axis){ //focntion inverse de index qui donne i ou j selon axe choisi
    int result;
    if (axis==0) 
    {
        result= I%Nx+1;
    }
    else if (axis==1)
    {
        result= I/Ny+1;
    }
    return result;
}

float maillage(int I,int axis){ // fonction qui donne xi ou xj selon axe choisi
    float result;
    if (axis==0) 
    {
         result=(couple(I,0)*1.0)/Nx;
    }
    else if (axis==1)
    {
        result=(couple(I,1)*1.0)/Ny;
    }
    return result;
    
}