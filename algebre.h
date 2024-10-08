#ifndef ALGEBRE_H
#define ALGEBRE_H









double* produit_MV( double* vecteur) ;

double produitScalaire(double *A, double *B, int taille);

void copierTableau(double *source, double *destination, int taille); 

void sommeTableaux(double *A, double *B, double *C, int taille);

void differenceTableaux(double *A, double *B, double *C, int taille);

void scalaireMultiplieTableau(double scalaire, double *vecteur, double *resultat, int taille);

double norm(double *A, int taille) ;

double* implicit_diff(double *source);







int indexe(int i,int j); // focntion qui donne l'indice de l'element stocké en fonction de i et j
int couple(int I,int axis); //focntion inverse de index qui donne i ou j selon axe choisi
double maillage(int I,int axis); // fonction qui donne xi ou xj selon axe choisi
 
    


#endif // ALGEBRE_H