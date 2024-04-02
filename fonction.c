#include "fonction.h"
#include <stdlib.h> 
#include "parametre.h"


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

                    for (int I=0,I<Nx*Ny,I++){








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
float* alpha( float x,float y){

}
float* beta( float x,float y){
    
}
float* gamma( float x,float y){
    
}