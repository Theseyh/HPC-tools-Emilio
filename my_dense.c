#include "spmv.h"

int my_dense(const unsigned int n, const double mat[], double vec[], double result[])
{
    // Initialiser le vecteur résultat à 0
    for (unsigned int i = 0; i < n; i++) {
        result[i] = 0.0;
    }

    // Multiplication matrice-vecteur
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            result[i] += mat[i * n + j] * vec[j];
        }
    }

    return 0; // Indiquer que l'exécution s'est bien passée
}

