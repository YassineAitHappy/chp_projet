#ifndef ALGEBRE_H
#define ALGEBRE_H









float* produit_MV( float* vecteur) ;










int indexe(int i,int j); // focntion qui donne l'indice de l'element stocké en fonction de i et j
int couple(int I,int axis); //focntion inverse de index qui donne i ou j selon axe choisi
float maillage(int I,int axis); // fonction qui donne xi ou xj selon axe choisi
 
    


#endif // ALGEBRE_H
