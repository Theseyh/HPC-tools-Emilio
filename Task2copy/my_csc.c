#include "spmv.h"
#include <stdio.h>

typedef struct {
    int row;   // Pointeur de début de chaque colonne
    int col;   // Indice de la ligne
    double val; // Valeur non-nulle
} SparseMat;

int my_csc(const unsigned int n, SparseMat sparse[], double vec[], double result[])
{

    // Parcourir chaqprinue colonne
    for (unsigned int j = 0; j < n; j++) {
        // Récupérer la plage d'indices pour cette colonne
        int start = sparse[j].col;
        int end = sparse[j+1].col;

        // Itérer sur les éléments non-nuls de cette colonne
        for (int k = start; k < end; k++) {
            int i = sparse[k].row;      // Indice de ligne
            double val = sparse[k].val; // Valeur non-nulle

            // Accumuler dans le résultat pour la ligne correspondante
            result[i] += val * vec[j];

        }
    }

    return 0; 
}
