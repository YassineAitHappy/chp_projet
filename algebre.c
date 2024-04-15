#include "algebre.h"
#include "math.h"
#include "parametre.h"
#include <stdlib.h> 
#include "fonction.h"
#include <stdio.h>
extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern int Nx, Ny, cas;


double* produit_MV( double* vecteur) {

    double* resultat = (double*)malloc(Nx*Ny* sizeof(double)); 


    // Calcul du produit matrice-vecteur
    if (vecteur== NULL) {
        
        return NULL;
    }
    else{

            switch(space_scheme) {


                // case 1://Centré
                //     for(int j=1; j<=Ny;j++){
                //         if(j==1){
                //             for(int i=1;i<=Nx;i++){
                //                 if(i==1){
                //                     resultat[indexe(i,j)]=alpha(maillage(indexe(i,j)+1,0),maillage(indexe(i,j)+1,1))*vecteur[indexe(i,j)+1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-alpha(maillage(indexe(i,j)+Nx-1,0),maillage(indexe(i,j)+Nx-1,1))*vecteur[indexe(i,j)+Nx-1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+beta(maillage(indexe(i,j)+Nx,0),maillage(indexe(i,j)+Nx,1))*vecteur[indexe(i,j)+Nx];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-beta(maillage(indexe(i,j)+Nx*(Ny-1),0),maillage(indexe(i,j)+Nx*(Ny-1),1))*vecteur[indexe(i,j)+Nx*(Ny-1)];
                //                     resultat[indexe(i,i)]=dt_imp*resultat[indexe(i,i)]+vecteur[indexe(i,i)];
                //                 }
                //                 else if(i==Nx){
                //                     resultat[indexe(i,j)]=alpha(maillage(indexe(i,j)-Nx+1,0),maillage(indexe(i,j)-Nx+1,1))*vecteur[indexe(i,j)-Nx+1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-alpha(maillage(indexe(i,j)-1,0),maillage(indexe(i,j)-1,1))*vecteur[indexe(i,j)-1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+beta(maillage(indexe(i,j)+Nx,0),maillage(indexe(i,j)+Nx,1))*vecteur[indexe(i,j)+Nx];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-beta(maillage(indexe(i,j)+Nx*(Ny-1),0),maillage(indexe(i,j)+Nx*(Ny-1),1))*vecteur[indexe(i,j)+Nx*(Ny-1)];
                //                     resultat[indexe(i,i)]=dt_imp*resultat[indexe(i,i)]+vecteur[indexe(i,i)];
                //                 }
                //                 else{
                //                     resultat[indexe(i,j)]=-alpha(maillage(indexe(i,j)-1,0),maillage(indexe(i,j)-1,1))*vecteur[indexe(i,j)-1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+alpha(maillage(indexe(i,j)+1,0),maillage(indexe(i,j)+1,1))*vecteur[indexe(i,j)+1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+beta(maillage(indexe(i,j)+Nx,0),maillage(indexe(i,j)+Nx,1))*vecteur[indexe(i,j)+Nx];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-beta(maillage(indexe(i,j)+(Ny-1)*Nx,0),maillage(indexe(i,j)+(Ny-1)*Nx,1))*vecteur[indexe(i,j)+(Ny-1)*Nx];
                //                     if(i==j){
                //                         resultat[indexe(i,i)]=dt_imp*resultat[indexe(i,i)]+vecteur[indexe(i,i)];
                //                     }
                //                 }
                //             }
                //         }
                //         else if(j==Ny){
                //             for(int i=1;i<=Nx;i++){
                //                 if(i==1){
                //                     resultat[indexe(i,j)]=beta(maillage(indexe(i,j)-Nx*(Ny-1),0),maillage(indexe(i,j)-Nx*(Ny-1),1))*vecteur[indexe(i,j)-Nx*(Ny-1)];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-beta(maillage(indexe(i,j)-Nx,0),maillage(indexe(i,j)-Nx,1))*vecteur[indexe(i,j)-Nx];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+alpha(maillage(indexe(i,j)+1,0),maillage(indexe(i,j)+1,1))*vecteur[indexe(i,j)+1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-alpha(maillage(indexe(i,j)+Nx-1,0),maillage(indexe(i,j)+Nx-1,1))*vecteur[indexe(i,j)+Nx-1];
                //                     resultat[indexe(i,i)]=dt_imp*resultat[indexe(i,i)]+vecteur[indexe(i,i)];
                //                 }
                //                 else if(i==Nx){
                //                     resultat[indexe(i,j)]=beta(maillage(indexe(i,j)-Nx*(Ny-1),0),maillage(indexe(i,j)-Nx*(Ny-1),1))*vecteur[indexe(i,j)-Nx*(Ny-1)];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-beta(maillage(indexe(i,j)-Nx,0),maillage(indexe(i,j)-Nx,1))*vecteur[indexe(i,j)-Nx];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+alpha(maillage(indexe(i,j)-Nx+1,0),maillage(indexe(i,j)-Nx+1,1))*vecteur[indexe(i,j)-Nx+1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-alpha(maillage(indexe(i,j)-1,0),maillage(indexe(i,j)-1,1))*vecteur[indexe(i,j)-1];
                //                     resultat[indexe(i,i)]=dt_imp*resultat[indexe(i,i)]+vecteur[indexe(i,i)];
                //                 }
                //                 else{
                //                     resultat[indexe(i,j)]=beta(maillage(indexe(i,j)-Nx*(Ny-1),0),maillage(indexe(i,j)-Nx*(Ny-1),1))*vecteur[indexe(i,j)-Nx*(Ny-1)];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-beta(maillage(indexe(i,j)-Nx,0),maillage(indexe(i,j)-Nx,1))*vecteur[indexe(i,j)-Nx];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+alpha(maillage(indexe(i,j)+1,0),maillage(indexe(i,j)+1,1))*vecteur[indexe(i,j)+1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-alpha(maillage(indexe(i,j)-1,0),maillage(indexe(i,j)-1,1))*vecteur[indexe(i,j)-1];
                //                     if(i==j){
                //                         resultat[indexe(i,i)]=dt_imp*resultat[indexe(i,i)]+vecteur[indexe(i,i)];
                //                     }
                //                 }
                //             }
                //         }
                        
                //         else if(j!=1 && j!=Ny){
                //             for(int i=1;i<=Nx;i++){
                //                 if(i%Nx==1){
                //                     resultat[indexe(i,j)]=-beta(maillage(indexe(i,j)-Nx,0),maillage(indexe(i,j)-Nx,1))*vecteur[indexe(i,j)-Nx];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+alpha(maillage(indexe(i,j)+1,0),maillage(indexe(i,j)+1,1))*vecteur[indexe(i,j)+1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-alpha(maillage(indexe(i,j)+Nx-1,0),maillage(indexe(i,j)+Nx-1,1))*vecteur[indexe(i,j)+Nx-1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+beta(maillage(indexe(i,j)+Nx,0),maillage(indexe(i,j)+Nx,1))*vecteur[indexe(i,j)+Nx];
                            
                //                     resultat[indexe(i,i)]=dt_imp*resultat[indexe(i,i)]+vecteur[indexe(i,i)];
                                    
                //                 }
                //                 else if(i%Nx==0){
                //                     resultat[indexe(i,j)]=-beta(maillage(indexe(i,j)-Nx,0),maillage(indexe(i,j)-Nx,1))*vecteur[indexe(i,j)-Nx];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+alpha(maillage(indexe(i,j)-Nx+1,0),maillage(indexe(i,j)-Nx+1,1))*vecteur[indexe(i,j)-Nx+1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-alpha(maillage(indexe(i,j)-1,0),maillage(indexe(i,j)-1,1))*vecteur[indexe(i,j)-1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+beta(maillage(indexe(i,j)+Nx,0),maillage(indexe(i,j)+Nx,1))*vecteur[indexe(i,j)+Nx];
                //                     resultat[indexe(i,i)]=dt_imp*resultat[indexe(i,i)]+vecteur[indexe(i,i)];
                //                 }
                //                 else{
                //                     resultat[indexe(i,j)]=-beta(maillage(indexe(i,j)-Nx,0),maillage(indexe(i,j)-Nx,1))*vecteur[indexe(i,j)-Nx];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+alpha(maillage(indexe(i,j)+1,0),maillage(indexe(i,j)+1,1))*vecteur[indexe(i,j)+1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]-alpha(maillage(indexe(i,j)-1,0),maillage(indexe(i,j)-1,1))*vecteur[indexe(i,j)-1];
                //                     resultat[indexe(i,j)]=resultat[indexe(i,j)]+beta(maillage(indexe(i,j)+Nx,0),maillage(indexe(i,j)+Nx,1))*vecteur[indexe(i,j)+Nx];
                //                     if(i==j){
                //                         resultat[indexe(i,i)]=dt_imp*resultat[indexe(i,i)]+vecteur[indexe(i,i)];
                //                     }
                //                 }
                //             }
                //         }

                        
                //     double* resultat_copy = (double*)malloc(Nx*Ny* sizeof(double));
                //     copierTableau(resultat,resultat_copy,Nx*Ny);
                //     scalaireMultiplieTableau(dt_imp,resultat,resultat,Nx*Ny);
                //     sommeTableaux(resultat,resultat_copy,resultat,Nx*Ny);

                    
                    
                //     }
                    
                //     break; 
                case 1:
                    for (int I=0;I<Nx*Ny;I++){
                        if(I==0){
                            resultat[I]=alpha(maillage(I+1,0),maillage(I+1,1))*vecteur[I+1]-alpha(maillage(I+Nx-1,0),maillage(I+Nx-1,1))*vecteur[I+Nx-1]+beta(maillage(I+Nx,0),maillage(I+Nx,1))*vecteur[I+Nx]-beta(maillage(I+Nx*(Ny-1),0),maillage(I+Nx*(Ny-1),1))*vecteur[I+Nx*(Ny-1)];
                        }
                        else if((0<I)&&(I<Nx-1)){
                            resultat[I]=-alpha(maillage(I-1,0),maillage(I-1,1))*vecteur[I-1]+alpha(maillage(I+1,0),maillage(I+1,1))*vecteur[I+1]+beta(maillage(I+Nx,0),maillage(I+Nx,1))*vecteur[I+Nx]-beta(maillage(I+(Ny-1)*Nx,0),maillage(I+(Ny-1)*Nx,1))*vecteur[I+(Ny-1)*Nx];
                            
                        }
                        else if(I==Nx-1){
                            resultat[I]=alpha(maillage(I-Nx+1,0),maillage(I-Nx+1,1))*vecteur[I-Nx+1]-alpha(maillage(I-1,0),maillage(I-1,1))*vecteur[I-1]+beta(maillage(I+Nx,0),maillage(I+Nx,1))*vecteur[I+Nx]-beta(maillage(I+Nx*(Ny-1),0),maillage(I+Nx*(Ny-1),1))*vecteur[I+Nx*(Ny-1)];
                        }
                        else if(I==Nx*(Ny-1)){
                            resultat[I]=beta(maillage(I-Nx*(Ny-1),0),maillage(I-Nx*(Ny-1),1))*vecteur[I-Nx*(Ny-1)]-beta(maillage(I-Nx,0),maillage(I-Nx,1))*vecteur[I-Nx]+alpha(maillage(I+1,0),maillage(I+1,1))*vecteur[I+1]-alpha(maillage(I+Nx-1,0),maillage(I+Nx-1,1))*vecteur[I+Nx-1];
                        }
                        else if(I==Nx*(Ny-1)+Nx-1){
                            resultat[I]=beta(maillage(I-Nx*(Ny-1),0),maillage(I-Nx*(Ny-1),1))*vecteur[I-Nx*(Ny-1)]-beta(maillage(I-Nx,0),maillage(I-Nx,1))*vecteur[I-Nx]+alpha(maillage(I-Nx+1,0),maillage(I-Nx+1,1))*vecteur[I-Nx+1]-alpha(maillage(I-1,0),maillage(I-1,1))*vecteur[I-1];

                        }
                        else if((Nx*(Ny-1)<I)&&(I<Nx*(Ny-1)+Nx-1)){
                            resultat[I]=beta(maillage(I-Nx*(Ny-1),0),maillage(I-Nx*(Ny-1),1))*vecteur[I-Nx*(Ny-1)]-beta(maillage(I-Nx,0),maillage(I-Nx,1))*vecteur[I-Nx]+alpha(maillage(I+1,0),maillage(I+1,1))*vecteur[I+1]-alpha(maillage(I-1,0),maillage(I-1,1))*vecteur[I-1];
                            
                        }
                        else if((Nx-1<I)&&(I<Nx*(Ny-1))){
                            
                            if(I%Nx==0){
                                    resultat[I]= -beta(maillage(I-Nx,0),maillage(I-Nx,1))*vecteur[I-Nx]+alpha(maillage(I+1,0),maillage(I+1,1))*vecteur[I+1]-alpha(maillage(I+Nx-1,0),maillage(I+Nx-1,1))*vecteur[I+Nx-1]+beta(maillage(I+Nx,0),maillage(I+Nx,1))*vecteur[I+Nx];
                                }
                            else if(I%Nx==Nx-1){
                                    resultat[I]=-beta(maillage(I-Nx,0),maillage(I-Nx,1))*vecteur[I-Nx]+alpha(maillage(I-Nx+1,0),maillage(I-Nx+1,1))*vecteur[I-Nx+1]-alpha(maillage(I-1,0),maillage(I-1,1))*vecteur[I-1]+beta(maillage(I+Nx,0),maillage(I+Nx,1))*vecteur[I+Nx];
                                    
                                }
                            
                            else if((I%Nx!=0) && (I%Nx!=Nx-1)){
                                resultat[I]=-beta(maillage(I-Nx,0),maillage(I-Nx,1))*vecteur[I-Nx]+alpha(maillage(I+1,0),maillage(I+1,1))*vecteur[I+1]-alpha(maillage(I-1,0),maillage(I-1,1))*vecteur[I-1]+beta(maillage(I+Nx,0),maillage(I+Nx,1))*vecteur[I+Nx];
                                
                                    
                            }
                        }
                        
                        
                    }
                    

                    break;

                case 2://Upwind

                    for (int I=0;I<Nx*Ny;I++){
                        
                            
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
        result= I/Nx+1;
    }
    return result;
}

double maillage(int I,int axis){ // fonction qui donne xi ou xj selon axe choisi
    double result;
    if (axis==0) 
    {
         result=xmin+couple(I,0)*dx;
    }
    else if (axis==1)
    {
        result=ymin+couple(I,1)*dy;
    }
    return result;
    
}

double produitScalaire(double *A, double *B, int taille) {
    double produit = 0;
    for (int i = 0; i < taille; i++) {
        produit += A[i] * B[i];
    }
    return produit;
}
double norm(double *A, int taille) {
    double produit = 0.0;
    double resultat ;
    for (int i = 0; i < taille; i++) {
        produit += A[i] * A[i];
    }
    resultat=sqrt(produit);
    return resultat;
}

void copierTableau(double *source, double *destination, int taille) {
    for (int i = 0; i < taille; i++) {
        destination[i] = source[i];
    }
}

void sommeTableaux(double *A, double *B, double *C, int taille) {
    for (int i = 0; i < taille; i++) {
        C[i] = A[i] + B[i];
    }
}

void differenceTableaux(double *A, double *B, double *C, int taille) {
    for (int i = 0; i < taille; i++) {
        C[i] = A[i] - B[i];
    }
}

void scalaireMultiplieTableau(double scalaire, double *vecteur, double *resultat, int taille) {
    for (int i = 0; i < taille; i++) {
        resultat[i] = scalaire * vecteur[i];
    }
}

double* implicit_diff(double *source) { //effectue un+1-dt*A*un+1 utilisé en bicgstab


    double* resultat = (double*)malloc(Nx*Ny* sizeof(double)); 



                    copierTableau(produit_MV(source),resultat,Nx*Ny);
                    scalaireMultiplieTableau(dt_imp,resultat,resultat,Nx*Ny);
                    differenceTableaux(source,resultat,resultat,Nx*Ny);
                    return resultat;
}
                