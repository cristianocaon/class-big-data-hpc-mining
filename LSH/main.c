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

    int d = 0, ii = 0, jj = 0, kk = 0, first = 1, match = 0, nclusters = 1, count = 1;;

    double inner_product = 0.0;

    int* data_hashes = malloc(sizeof(int) * ndata * m);
    int* cluster_assign = malloc(sizeof(int) * ndata);

    for (int i = 0; i < ndata * dim; i += dim) {

        jj = ii;
        kk = ii;

        for (int j = 0; j < m; j++) {
            for (int k = i; k < i + dim; k++) {
                inner_product += data[k] * h[j][d];
                d += 1;
            }
            data_hashes[ii] = (int)floor((inner_product - b[j]) / W);
            inner_product = 0.0;
            ii += 1;
            d = 0;
        }

        if (first) {
            cluster_assign[0] = 0;
            cluster_hashval[0] = malloc(sizeof(int) * m);
            for (int j = 0; j < m; j++) {
                cluster_hashval[0][j] = data_hashes[kk];
                kk += 1;
            }
            first = 0;
            continue;
        }

        kk = jj;

        for (int j = 0; j < nclusters; j++) {
            for (int k = 0; k < m; k++) {
                if (data_hashes[kk] == cluster_hashval[j][k]) {
                    match += 1;
                    kk += 1;
                }
                else {
                    match = 0;
                    kk -= k;
                    break;
                }
            }

            // Assign data point to cluster with same hash
            if (match == m) {
                printf("\nCluster matched!");
                cluster_assign[count] = j;
                break;
            }

            kk = jj;
        }

        // No cluster found for this hash, create new cluster.
        if (match == 0) {
            printf("\nNo cluster matched, creating new one...");
            cluster_hashval = realloc(cluster_hashval, sizeof(*cluster_hashval) * (nclusters + 1));
            cluster_hashval[nclusters] = malloc(sizeof(int) * m);
            cluster_assign[count] = nclusters;
            for (int j = 0; j < m; j++) {
                cluster_hashval[nclusters][j] = data_hashes[kk];
                kk += 1;
            }
            nclusters += 1;
        }

        count += 1;
    }

    printf("\n\n");
    for (int i = 0; i < ndata * m; i++) {
        printf("%d\t", data_hashes[i]);
        if ((i + 1) % m == 0) {
            printf("\n");
        }
    }

    for (int i = 0; i < nclusters; i++) {
        printf("\nCluster %d hashes: ", i);
        for (int j = 0; j < m; j++) {
            printf("%d\t", cluster_hashval[i][j]);
        }
    }

    for (int i = 0; i < ndata; i++) {
        printf("\n%d", cluster_assign[i]);
    }

    free(data_hashes);
    free(cluster_assign);

    return nclusters;
}


int search_LSH(int dim, int ndata, double* data,
    int m, double W, double** h, double* b,
    int nclusters, int* cluster_start, int* cluster_size, int** cluster_hashval,
    double* query_pt, double* result_pt) {
    /*******************************************************************
    * The function returns...
    *
    ******************************************************************/

    int ii = 0, checked = 0, match = 0, cluster_id = -1;
    int* query_hash = malloc(sizeof(int) * m);

    double inner_product = 0.0, dist = 0.0, min_dist = __DBL_MAX__;
    double* closest = malloc(sizeof(double) * dim);

    // Computing query point hash.
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < dim; j++) {
            inner_product += query_pt[j] * h[i][j];
        }
        query_hash[i] = (int)floor((inner_product - b[i]) / W);
        inner_product = 0.0;
    }

    // Searching cluster for query point.
    for (int i = 0; i < nclusters; i++) {
        for (int j = 0; j < m; j++) {
            if (query_hash[j] == cluster_hashval[i][j]) {
                match += 1;
            }
            else {
                match = 0;
                break;
            }
        }

        // Found cluster for query point.
        if (match == m) {
            printf("\nFound cluster for query point!\n");
            cluster_id = i;
            break;
        }
    }

    // Did not find cluster for query point.
    if (match == 0) {
        printf("\nDid not find cluster for query point.\n");
        exit(1);
    }

    // Looking for closest data point to query point inside cluster.
    for (int i = cluster_start[cluster_id]; i < cluster_start[cluster_id] + cluster_size[cluster_id] * dim; i += dim) {
        for (int j = i; j < i + dim; j++) {
            dist += fabs(query_pt[j]) - fabs(data[j]);
        }

        if (dist < min_dist) {
            min_dist = dist;
            for (int k = i; k < i + dim; k++) {
                closest[ii] = data[k];
                ii += 1;
            }
            ii = 0;
        }

        dist = 0.0;
        checked += 1;
    }

    printf("\nChecked %d data points for query point.", checked);

    printf("\nCloset data point: ");
    for (int i = 0; i < dim; i++) {
        printf("%f\t", closest[i]);
    }

    free(closest);
    free(query_hash);

    return checked;
}


int main() {
    /********************************************************************************
     * Steps performed in the program:
     *
    ********************************************************************************/

    // int ndata = 1000000, dim = 16, m = 5, nclusters = 1;
    int ndata = 10, dim = 4, m = 5, nclusters = 1;

    int* cluster_size = malloc(sizeof(int) * nclusters);
    int* cluster_start = malloc(sizeof(int) * nclusters);
    int** cluster_hashval = malloc(sizeof(*cluster_hashval) * nclusters * m);

    double inner_product, sum = 0.0, W = 0.3;

    double* centroids = malloc(sizeof(double) * dim);
    double* data = malloc(sizeof(double) * ndata * dim);
    double* b = malloc(sizeof(double) * m);
    double** h = malloc(sizeof(double*) * m * dim);

    printf("\nParameters:\n+----------------------+");
    printf("\nndata = %d\ndim = %d\nm = %d\nW = %.2f\n", ndata, dim, m, W);
    printf("+----------------------+\n\n");

    printf("\nGenerating random data...\n");
    for (int i = 0; i < ndata * dim; i++) {
        data[i] = (double)rand() / RAND_MAX;
    }

    printf("\nGenerating hashes...\n");
    for (int i = 0; i < m; i++) {
        h[i] = malloc(sizeof(double) * dim);
        for (int j = 0; j < dim; j++) {
            h[i][j] = 2.0 * (rand() / RAND_MAX) - 1.0;
        }
    }

    printf("\nCalculating centroids...\n");
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < ndata * dim; j += dim) {
            sum += data[j];
        }
        centroids[i] = sum / ndata;
        sum = 0.0;
    }

    printf("\nGenerating b values...\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < dim; j++) {
            inner_product += centroids[i] * h[i][j];
        }
        b[i] = inner_product;
        inner_product = 0.0;
    }

    printf("\nSorting data using LSH...\n");
    LSH(dim, ndata, data, m, W, h, b, cluster_start, cluster_size, cluster_hashval);

    free(b);
    free(h);
    free(data);
    free(centroids);
    free(cluster_size);
    free(cluster_start);
    // free(cluster_hashval);

    return 0;
}