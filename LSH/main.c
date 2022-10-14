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
    int** cluster_start, int** cluster_size, int*** cluster_hashval) {
    /*******************************************************************
     *  Organizes the data array using the LSH algorithm.
     *  Returns number of clusters created.
     *
     *  Array sizes: h[m][dim], b[m], cluster_hashval[nclusters][m]
     *      b[i] may choose <x_c, h_i>, where x_c [dim] is centroid.
     *
     *  Note: You don't have value of nclusters as input. It's output.
     *
     *  Question:
     *      How to keep cost within O(nclusters * m * ndata * dim)?
     *      Your code must have this computation complexity.
    ******************************************************************/

    int first = 1, match = 0, count = 1, nclusters = 1;
    int* data_hash = (int*)malloc(sizeof(int) * m);

    double inner_product = 0.0;
    double** clustered_data = (double**)malloc(sizeof(double*));

    // Iterating through all data points.
    for (int i = 0; i < ndata * dim; i += dim) {
        // Calculating data point hash.
        for (int j = 0; j < m; j++) {
            for (int k = i, ii = 0; k < i + dim; k++, ii++) {
                inner_product += data[k] * h[j][ii];
            }
            data_hash[j] = (int)floor((inner_product - b[j]) / W);
            inner_product = 0.0;
        }

        // Initializing variables for first data point only.
        if (first) {
            (*cluster_size)[0] = 1;
            (*cluster_hashval)[0] = (int*)malloc(sizeof(int) * m);
            for (int j = 0; j < m; j++) {
                (*cluster_hashval)[0][j] = data_hash[j];
            }
            clustered_data[0] = (double*)malloc(sizeof(double) * dim);
            for (int j = 0; j < dim; j++) {
                clustered_data[0][j] = data[j];
            }
            first = 0;
            continue;
        }

        // Comparing data point hash to existing cluster hashes.
        for (int j = 0; j < nclusters; j++) {
            for (int k = 0; k < m; k++) {
                if (data_hash[k] == (*cluster_hashval)[j][k]) {
                    match += 1;
                }
                else {
                    match = 0;
                    break;
                }
            }

            // Data point hash matched cluster hash. 
            // Assign data point to cluster and update variables.
            if (match == m) {
                clustered_data[j] = (double*)realloc(clustered_data[j], sizeof(double) * dim * ((*cluster_size)[j] + 1));
                for (int k = i, ii = (*cluster_size)[j] * dim; k < i + dim; k++, ii++) {
                    clustered_data[j][ii] = data[k];
                }
                (*cluster_size)[j] += 1;
                break;
            }
        }

        // Data point hash did not match existing clusters' hashes. 
        // Create new cluster with data point hash and update variables.
        if (match == 0) {
            *cluster_hashval = (int**)realloc(*cluster_hashval, sizeof(int*) * (nclusters + 1));
            (*cluster_hashval)[nclusters] = (int*)malloc(sizeof(int) * m);
            for (int j = 0; j < m; j++) {
                (*cluster_hashval)[nclusters][j] = data_hash[j];
            }
            *cluster_size = (int*)realloc(*cluster_size, sizeof(int) * (nclusters + 1));
            (*cluster_size)[nclusters] = 1;
            clustered_data = (double**)realloc(clustered_data, sizeof(double*) * (nclusters + 1));
            clustered_data[nclusters] = (double*)malloc(sizeof(double) * dim);
            for (int j = i, k = 0; j < i + dim; j++, k++) {
                clustered_data[nclusters][k] = data[j];
            }
            nclusters += 1;
        }
        count += 1;
    }

    // Acquiring clusters start variables.
    *cluster_start = (int*)realloc(*cluster_start, sizeof(int) * nclusters);
    (*cluster_start)[0] = 0;
    for (int i = 1; i < nclusters; i++) {
        (*cluster_start)[i] = (*cluster_start)[i - 1] + (*cluster_size)[i - 1] * dim;
    }

    // Organizing data array in order of clusters.
    for (int i = 0, k = 0; i < nclusters; i++) {
        for (int j = 0; j < (*cluster_size)[i] * dim; j++, k++) {
            data[k] = clustered_data[i][j];
        }
    }

    // Clean up.
    for (int i = 0; i < nclusters; i++) {
        free(clustered_data[i]);
    }
    free(data_hash);
    free(clustered_data);

    return nclusters;
}


int search_LSH(int dim, int ndata, double* data,
    int m, double W, double** h, double* b,
    int nclusters, int* cluster_start, int* cluster_size, int** cluster_hashval,
    double* query_pt, double** result_pt) {
    /***************************************************************************
     * Searches closest data point to query point in LSH-organized data array.
     * Returns number of data points checked in data array.
    ***************************************************************************/

    int checked = 0, match = 0, cluster_id = -1;
    int* query_hash = malloc(sizeof(int) * m);

    double inner_product = 0.0, dist = 0.0, min_dist = __DBL_MAX__;

    // Computing query point hash.
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < dim; j++) {
            inner_product += query_pt[j] * h[i][j];
        }
        query_hash[i] = (int)floor((inner_product - b[i]) / W);
        inner_product = 0.0;
    }

    // Searching if query point belongs to a cluster.
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

        // Found query point inside a cluster.
        if (match == m) {
            printf("\nFound query point inside cluster %d!", i);
            cluster_id = i;
            break;
        }
    }

    // Query point is not inside any cluster.
    if (match == 0) {
        printf("\nQuery point is not inside any cluster. Returning...");
        return -1;
    }

    // Searching for closest data point to query point inside cluster.
    for (int i = cluster_start[cluster_id]; i < cluster_start[cluster_id] + cluster_size[cluster_id] * dim; i += dim) {
        for (int j = i, k = 0; j < i + dim; j++, k++) {
            dist += fabs(query_pt[k]) - fabs(data[j]);
        }
        // Found new closest data point.
        if (dist < min_dist) {
            min_dist = dist;
            for (int j = i, k = 0; j < i + dim; j++, k++) {
                (*result_pt)[k] = data[j];
            }
        }
        dist = 0.0;
        checked += 1;
    }

    printf("\nClosest data point from cluster %d: ", cluster_id);
    for (int i = 0; i < dim; i++) {
        printf("%f\t", (*result_pt)[i]);
    }

    // Clean up.
    free(query_hash);

    return checked;
}


int main() {
    /********************************************************************************
     * Steps performed in the program:
     *  1- Generate necessary parameters.
     *  2- Organize data array with given parameters using LSH.
     *  3- Search for the closest data point to query point in organized data array.
    ********************************************************************************/

    int nclusters, ndata = 1000000, dim = 16, m = 5, checked = 0;
    double inner_product = 0.0, sum = 0.0, W = 0.3;

    int* cluster_size = (int*)malloc(sizeof(int));
    int* cluster_start = (int*)malloc(sizeof(int));
    int** cluster_hashval = (int**)malloc(sizeof(int*));

    double* centroids = (double*)malloc(sizeof(double) * dim);
    double* data = (double*)malloc(sizeof(double) * ndata * dim);
    double* b = (double*)malloc(sizeof(double) * m);
    double* query_pt = (double*)malloc(sizeof(double) * dim);
    double* result_pt = (double*)malloc(sizeof(double) * dim);
    double** h = (double**)malloc(sizeof(double*) * m);

    printf("\nParameters:\n==============================");
    printf("\nndata = %d\ndim = %d\nm = %d\nW = %.2f\n", ndata, dim, m, W);
    printf("==============================\n");

    printf("\nGenerating random data...");
    for (int i = 0; i < ndata * dim; i++) {
        data[i] = (double)rand() / RAND_MAX;
    }

    printf("\nGenerating random h values...");
    for (int i = 0; i < m; i++) {
        h[i] = (double*)malloc(sizeof(double) * dim);
        for (int j = 0; j < dim; j++) {
            h[i][j] = (double)rand() / RAND_MAX * 2.0 - 1.0;
        }
    }

    printf("\nCalculating centroids from data dimensions...");
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < ndata * dim; j += dim) {
            sum += data[j];
        }
        centroids[i] = sum / ndata;
        sum = 0.0;
    }

    printf("\nGenerating b values with h and centroids...");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < dim; j++) {
            inner_product += centroids[i] * h[i][j];
        }
        b[i] = inner_product;
        inner_product = 0.0;
    }

    printf("\n\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n");
    printf("\nOrganizing data using LSH...");

    nclusters = LSH(dim, ndata, data, m, W, h, b, &cluster_start, &cluster_size, &cluster_hashval);

    printf("\nPartitioned data into %d clusters!", nclusters);

    printf("\n\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n");

    printf("\nGenerating query point...");
    for (int i = 0; i < dim; i++) {
        query_pt[i] = (double)rand() / RAND_MAX;
    }

    printf("\nSearching closest data point to query point...");
    checked = search_LSH(dim, ndata, data, m, W, h, b, nclusters, cluster_start, cluster_size, cluster_hashval, query_pt, &result_pt);

    if (checked != -1) {
        printf("\nChecked %d data points for query point.", checked);
    }

    printf("\nCleaning up...");

    for (int i = 0; i < m; i++) {
        free(h[i]);
    }
    for (int i = 0; i < nclusters; i++) {
        free(cluster_hashval[i]);
    }
    free(b);
    free(h);
    free(data);
    free(query_pt);
    free(result_pt);
    free(centroids);
    free(cluster_size);
    free(cluster_start);
    free(cluster_hashval);

    return 0;
}
