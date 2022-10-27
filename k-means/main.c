/******************************************

 * Texas Tech University
 * CS 5331-001: Mining Big Data with HPC
 * Professor: Yu Zhuang
 * Assignment 3: K-Means
 * Student: Cristiano Eleutherio Caon
 * R#: 11474435

******************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/******************************************************************************
 k :                        number of clusters, i.e. the K in K-mean.
 cluster_centroid[k*dim]:   input  -- stores initial k centers
                            output -- stores final k centers
 cluster_radius[k]:         output
                            Stores the radius of each output cluster
 cluster_start[k]:          output
                            Stores the index of each cluster's starting datum
 cluster_size[k]:           output
                            Stores the num of data points in each cluster
 cluster_assign[ndata]:     buffer that stores the membership of each datum
*******************************************************************************/


int initial_centers(int dim, int ndata, double* data, int kk, double*** cluster_centroid) {
    /********************************************************************************
     *
    ********************************************************************************/

    /*
    choose the 1st center randomly
    choose the datum furthest away from center1 as center2
    J=2
    while (J < K)
        for 𝑑𝑎𝑡𝑢𝑚_i looping through all data
            calculate the distance of 𝑑𝑎𝑡𝑢𝑚_i to all J centers
            the closest distance is stored into d_min[i]
        end of for-loop
        choose the datum whose d_min[i*] is the largest among all i from 1 to n
        datum[i*] is chosen as the (J+1)-th center
        increment J
    end of while-loop
    */

    int largest_idx = 0, clusters = 2;
    double temp_dist, largest, dist = __DBL_MIN__;

    double* d_min = (double*)malloc(sizeof(double) * ndata);

    // Generating first cluster centroid.
    (*cluster_centroid)[0] = (double*)malloc(sizeof(double) * dim);
    for (int i = 0; i < dim; i++) {
        (*cluster_centroid)[0][i] = data[i];
    }

    // Finding second furthest away centroid.
    (*cluster_centroid)[1] = (double*)malloc(sizeof(double) * dim);
    for (int i = 0; i < ndata * dim; i += dim) {
        temp_dist = 0.0;
        for (int j = i, k = 0; j < i + dim; j++, k++) {
            temp_dist += pow(fabs((*cluster_centroid)[0][k]) - fabs(data[j]), 2);
        }
        temp_dist = sqrt(temp_dist);
        if (temp_dist > dist) {
            for (int j = i, k = 0; j < i + dim; j++, k++) {
                (*cluster_centroid)[1][k] = data[j];
            }
            dist = temp_dist;
        }
    }

    // Creating new cluster centers.
    while (clusters < kk) {
        // Getting each data point closest center distance.
        for (int i = 0, n = 0; i < ndata * dim; i += dim, n++) {
            d_min[n] = __DBL_MAX__;
            for (int j = 0; j < clusters; j++) {
                temp_dist = 0.0;
                for (int k = i, ii = 0; k < i + dim; k++, ii++) {
                    temp_dist += pow(fabs((*cluster_centroid)[j][ii]) - fabs(data[k]), 2);
                }
                temp_dist = sqrt(temp_dist);
                if (temp_dist < d_min[n]) {
                    d_min[n] = temp_dist;
                }
            }
        }
        // Getting largest data point to a center distance.
        largest = d_min[0];
        for (int i = 1; i < ndata; i++) {
            if (d_min[i] > largest) {
                largest = d_min[i];
                largest_idx = i;
            }
        }
        // Saving data point as a new cluster center.
        (*cluster_centroid)[clusters] = (double*)malloc(sizeof(double) * dim);
        for (int i = largest_idx * dim, j = 0; i < largest_idx * dim + dim; i++, j++) {
            (*cluster_centroid)[clusters][j] = data[i];
        }
        clusters++;
    }

    return 0;
}


double kmeans(int dim, int ndata, double** data, int kk,   // input
    short** cluster_assign,                                // buffer
    int** cluster_start, int** cluster_size,               // output
    double** cluster_radius, double*** cluster_centroid) { // output
    /********************************************************************************
     * Returns the sum of square errors.
    ********************************************************************************/

    /*
     Choose K initial cluster centers
     stop_iteration = 0
     while (stop_iteration == 0) {
        count_cluster_change = 0
        for 𝑑𝑎𝑡𝑢𝑚_i looping through all data
            calculate the distance from 𝑑𝑎𝑡𝑢𝑚_i to all K cluster centers
            assign 𝑑𝑎𝑡𝑢𝑚_i to the clusters with the closest cluster center
            if (𝑑𝑎𝑡𝑢𝑚_i changes cluster assignment)
                increment count_cluster_change
        end of for-loop // comment: some data are assigned to different clusters
        calculate the centroid of each cluster, the new K centroid are new centers
        if (count_cluster_change == 0)
            stop_iteration = 1
     end of while-loop
    */

    int count, count_cluster_change, chosen_cluster, stop_iteration = 0;

    double temp_dist, d_min;

    double* buffer;

    while (stop_iteration == 0) {
        count_cluster_change = 0;
        for (int i = 0, n = 0; i < ndata * dim; i += dim, n++) {
            d_min = __DBL_MAX__;
            chosen_cluster = (*cluster_assign)[n];
            // Finding closest cluster center to data point.
            for (int j = 0; j < kk; j++) {
                temp_dist = 0.0;
                for (int k = i, ii = 0; k < i + dim; k++, ii++) {
                    temp_dist += pow(fabs((*cluster_centroid)[j][ii]) - fabs((*data)[k]), 2);
                }
                temp_dist = sqrt(temp_dist);
                if (temp_dist < d_min) {
                    d_min = temp_dist;
                    chosen_cluster = j;
                }
            }
            // Data point changed cluster assignment.
            if (chosen_cluster != (*cluster_assign)[n]) {
                (*cluster_assign)[n] = chosen_cluster;
                count_cluster_change++;
            }
        }

        // Re-calculating cluster centers.
        for (int i = 0; i < kk; i++) {
            for (int j = 0, k = 0; j < ndata * dim; j += dim, k++) {
                if ((*cluster_assign)[k] == i) {
                    printf("\nRecalculate center");
                }
            }
        }

        // Exit condition.
        stop_iteration = (count_cluster_change == 0) ? 1 : 0;
    }

    count = 0;
    (*cluster_start)[0] = 0;
    buffer = (double*)malloc(sizeof(double) * ndata * dim);
    for (int i = 0; i < kk; i++) {
        for (int j = 0, k = 0; j < ndata * dim; j += dim, k++) {
            if ((*cluster_assign)[k] == i) {
                for (int ii = j; ii < j + dim; ii++) {
                    buffer[count] = (*data)[ii];
                    count++;
                }
                (*cluster_size)[i]++;
            }
        }
        if (i > 0) {
            (*cluster_start)[i] = (*cluster_start)[i - 1] + (*cluster_size)[i - 1] * dim;
        }
    }

    (*data) = buffer;

    // printf("\n\n");
    // for (int i = 0; i < kk; i++) {
    //     printf("Cluster start: %d\tCluster size: %d\n", (*cluster_start)[i], (*cluster_size)[i]);
    // }

    // printf("\n\n");
    // for (int i = 0; i < ndata; i++) {
    //     printf("%d\n", (*cluster_assign)[i]);
    // }

    // printf("\n\n");
    // for (int i = 0; i < ndata * dim; i++) {
    //     printf("%f\t", buffer[i]);
    //     if ((i + 1) % dim == 0) {
    //         printf("\n");
    //     }
    // }

    return 0;
}


int search_kmeans(int dim, int ndata, double* data, int kk,
    int* cluster_start, int* cluster_size,
    double* cluster_radius, double** cluster_centroid,
    double* query_pt, double* result_pt) {
    /********************************************************************************
     * Returns the number of data points checked.
    ********************************************************************************/

    return 0;
}


int main() {
    /********************************************************************************
     * Steps performed in main program:
     *
    ********************************************************************************/

    int ndata = 10, dim = 3, kk = (int)sqrt(ndata);

    int* cluster_start = (int*)malloc(sizeof(int) * kk);
    int* cluster_size = (int*)malloc(sizeof(int) * kk);

    short* cluster_assign = (short*)malloc(sizeof(short) * ndata);

    double* data = (double*)malloc(sizeof(double) * ndata * dim);
    double* cluster_radius = (double*)malloc(sizeof(double) * kk);
    double** cluster_centroid = (double**)malloc(sizeof(double*) * kk);

    printf("\nParameters:\n+----------------------+");
    printf("\nndata = %d\ndim = %d\nkk = %d\n", ndata, dim, kk);
    printf("+----------------------+\n\n");

    printf("\nGenerating random data...\n");
    for (int i = 0; i < ndata * dim; i++) {
        data[i] = (double)rand() / RAND_MAX;
        printf("%f\t", data[i]);
        if ((i + 1) % dim == 0) {
            printf("\n");
        }
    }

    for (int i = 0; i < ndata; i++) {
        cluster_assign[i] = -1;
    }

    printf("\nCreating initial cluster centers...\n");
    initial_centers(dim, ndata, data, kk, &cluster_centroid);

    printf("\nSorting data with k-means...\n");
    kmeans(dim, ndata, &data, kk, &cluster_assign, &cluster_start, &cluster_size, &cluster_radius, &cluster_centroid);

    // for (int i = 0; i < ndata * dim; i++) {
    //     printf("%f\t", data[i]);
    //     if ((i + 1) % dim == 0) {
    //         printf("\n");
    //     }
    // }

    // printf("\nAll centroids:\n");
    // for (int i = 0; i < kk; i++) {
    //     for (int j = 0; j < dim; j++) {
    //         printf("%f\t", cluster_centroid[i][j]);
    //     }
    //     printf("\n");
    // }

    // Clean up
    for (int i = 0; i < kk; i++) {
        free(cluster_centroid[i]);
    }
    free(data);
    free(cluster_size);
    free(cluster_start);
    free(cluster_assign);
    free(cluster_radius);

    return 0;
}
