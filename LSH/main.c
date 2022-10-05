/******************************************

 * Texas Tech University
 * CS 5331-001: Mining Big Data with HPC
 * Professor: Yu Zhuang
 * Assignment 2: LSH
 * Student: Cristiano Eleutherio Caon
 * R#: 11474435

******************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int LSH(int dim, int ndata, double* data,
    int m, double W, double** h, double* b,
    int* cluster_start, int* cluster_size, int** cluster_hashval) {
    /*******************************************************************
    * The function returns the number of buckets or clusters
    *
    * Array sizes: h[m][dim], b[m], cluster_hashval[nclusters][m]
    *           b[i] may choose <x_c, h_i>, where x_c [dim] is centroid
    *
    * Note: You don't have value of nclusters as input. It's output.
    *
    * Question:
    *      How to keep cost of within O(nclusters * m * ndata * dim)?
    *      Your code must have this computation complexity.
    ******************************************************************/

}


int search_LSH(int dim, int ndata, double* data,
    int m, double W, double** h, double* b,
    int nclusters, int* cluster_start, int* cluster_size, int** cluster_hashval,
    double* query_pt, double* result_pt) {
    /*******************************************************************
    * The function returns...
    *
    ******************************************************************/
}


int main() {
    /********************************************************************************
     * Steps performed in the program:
     *
    ********************************************************************************/

    int ndata = 1000000, dim = 16, m = (int)(dim / 2), W = 1;

    int* cluster_size = malloc(sizeof(int) * kk);
    int* cluster_start = malloc(sizeof(int) * kk);
    int* cluster_hashval = malloc(sizeof(int) * kk);

    double* b = malloc(sizeof(double) * m);
    double* h = malloc(sizeof(double) * m * dim);
    double* data = malloc(sizeof(double) * ndata * dim);

    printf("\nParameters:\n+----------------------+");
    printf("\nndata = %d\ndim = %d\nm = %d\nW = %d\n", ndata, dim, m, W);
    printf("+----------------------+\n\n");

    printf("\nGenerating random data...\n");
    for (int i = 0; i < ndata * dim; i++) {
        data[i] = (double)rand() / RAND_MAX;
    }

    printf("\nGenerating random hashes...\n");
    for (int i = 0; i < m * dim; i++) {
        h[i] = (double)rand() / RAND_MAX;
    }

    printf("\nGenerating random b values...\n");
    for (int i = 0; i < m; i++) {
        b[i] = (double)rand() / RAND_MAX;
    }

    printf("\nSorting data using LSH...\n");

    free(data);
    free(cluster_size);
    free(cluster_start);
    free(cluster_hashval);

    return 0;
}