float produitScalaire(float *A, float *B, int taille) {
    float produit = 0;
    for (int i = 0; i < taille; i++) {
        produit += A[i] * B[i];
    }
    return produit;
}

