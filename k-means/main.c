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


int initial_centers(int dim, int ndata, double* data, int kk, double*** cluster_centroid) {
    /******************************************************
     * Returns the initial centers for the 'kk' clusters.
    ******************************************************/

    int largest_idx = 0, clusters = 2;
    double temp_dist, largest, dist = __DBL_MIN__;

    double* d_min = (double*)malloc(sizeof(double) * ndata);

    // Generating first cluster centroid
    (*cluster_centroid)[0] = (double*)malloc(sizeof(double) * dim);
    for (int i = 0; i < dim; i++) {
        (*cluster_centroid)[0][i] = data[i];
    }

    // Finding second furthest away centroid
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

    // Creating new cluster centers
    while (clusters < kk) {
        // Getting each data point closest center distance
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
        // Getting largest data point to a center distance
        largest = d_min[0];
        for (int i = 1; i < ndata; i++) {
            if (d_min[i] > largest) {
                largest = d_min[i];
                largest_idx = i;
            }
        }
        // Saving data point as a new cluster center
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
    /****************************************************
     * Organizes the data array using K-Means algorithm
     *  into 'kk' number of clusters.
     *
     * Returns the sum of square errors.
    ****************************************************/

    int count, count_cluster_change, chosen_cluster, stop_iteration = 0;

    double temp_dist, d_min, total_error, cluster_error, point_error;

    double* data_buffer;
    double* centroid_buffer;

    while (stop_iteration == 0) {
        count_cluster_change = 0;
        for (int i = 0, n = 0; i < ndata * dim; i += dim, n++) {
            d_min = __DBL_MAX__;
            chosen_cluster = (*cluster_assign)[n];
            // Finding closest cluster center to data point
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
            // Data point changed cluster assignment
            if (chosen_cluster != (*cluster_assign)[n]) {
                (*cluster_assign)[n] = chosen_cluster;
                count_cluster_change++;
            }
        }

        // Exit condition
        stop_iteration = (count_cluster_change == 0) ? 1 : 0;

        // Re-calculating cluster centers
        for (int i = 0; i < kk; i++) {
            count = 0;
            centroid_buffer = (double*)malloc(sizeof(double) * dim);
            for (int j = 0, k = 0; j < ndata * dim; j += dim, k++) {
                if ((*cluster_assign)[k] == i) {
                    for (int ii = 0; ii < dim; ii++) {
                        centroid_buffer[ii] += (*data)[j + ii];
                    }
                    count++;
                }
            }
            for (int j = 0; j < dim; j++) {
                centroid_buffer[j] /= count;
                (*cluster_centroid)[i][j] = centroid_buffer[j];
            }
        }

        // Re-calculating cluster radius
        for (int i = 0; i < kk; i++) {
            (*cluster_radius)[i] = __DBL_MIN__;
            for (int j = 0, k = 0; j < ndata * dim; j += dim, k++) {
                if ((*cluster_assign)[k] == i) {
                    temp_dist = 0.0;
                    for (int ii = 0; ii < dim; ii++) {
                        temp_dist += pow(fabs((*cluster_centroid)[i][ii]) - fabs((*data)[j + ii]), 2);
                    }
                    temp_dist = sqrt(temp_dist);
                    (*cluster_radius)[i] = (temp_dist > (*cluster_radius)[i]) ? temp_dist : (*cluster_radius)[i];
                }
            }
        }
    }

    // Sorting data array in order of cluster data points
    count = 0;
    (*cluster_start)[0] = 0;
    data_buffer = (double*)malloc(sizeof(double) * ndata * dim);
    for (int i = 0; i < kk; i++) {
        for (int j = 0, k = 0; j < ndata * dim; j += dim, k++) {
            if ((*cluster_assign)[k] == i) {
                for (int ii = j; ii < j + dim; ii++) {
                    data_buffer[count] = (*data)[ii];
                    count++;
                }
                (*cluster_size)[i]++;
            }
        }
        if (i > 0) {
            (*cluster_start)[i] = (*cluster_start)[i - 1] + (*cluster_size)[i - 1] * dim;
        }
    }
    (*data) = data_buffer;

    // Computing the sum of square errors
    total_error = 0.0;
    for (int i = 0; i < kk; i++) {
        cluster_error = 0.0;
        for (int j = (*cluster_start)[i]; j < (*cluster_start)[i] + (*cluster_size)[i] * dim; j += dim) {
            point_error = 0.0;
            for (int k = j, ii = 0; k < dim; k++, ii++) {
                point_error += fabs((*cluster_centroid)[k][ii]) - fabs((*data)[k]);
            }
            cluster_error += pow(point_error, 2);
        }
        total_error += cluster_error / (*cluster_size)[i];
    }

    // printf("\n\n");
    // for (int i = 0; i < kk; i++) {
    //     printf("Cluster start: %d\tCluster size: %d\tCluster radius: %f\n", (*cluster_start)[i], (*cluster_size)[i], (*cluster_radius)[i]);
    // }

    // printf("\n\n");
    // for (int i = 0; i < ndata; i++) {
    //     printf("%d\n", (*cluster_assign)[i]);
    // }

    // printf("\n\n");
    // for (int i = 0; i < ndata * dim; i++) {
    //     printf("%f\t", data_buffer[i]);
    //     if ((i + 1) % dim == 0) {
    //         printf("\n");
    //     }
    // }

    free(centroid_buffer);

    return total_error;
}


int search_kmeans(int dim, int ndata, double* data, int kk,
    int* cluster_start, int* cluster_size,
    double* cluster_radius, double** cluster_centroid,
    double* query_pt, double* result_pt) {
    /*****************************************************
     * Searches for closest data point from data array to
     *  given query points.
     *
     * Returns the number of data points checked
    *****************************************************/

    return 0;
}


int main() {
    /***********************************************
     * Steps performed in main program:
     *
     *  1- Generate random data
     *  2- Create initial cluster centers
     *  3- Sort data with K-Means
     *  4- Search closest distance to query points
    ***********************************************/

    int ndata = 10000, dim = 16, kk = (int)sqrt(ndata);

    int* cluster_start = (int*)malloc(sizeof(int) * kk);
    int* cluster_size = (int*)malloc(sizeof(int) * kk);

    short* cluster_assign = (short*)malloc(sizeof(short) * ndata);

    double error;

    double* data = (double*)malloc(sizeof(double) * ndata * dim);
    double* cluster_radius = (double*)malloc(sizeof(double) * kk);
    double** cluster_centroid = (double**)malloc(sizeof(double*) * kk);

    printf("\nParameters:\n+----------------------+");
    printf("\nndata = %d\ndim = %d\nkk = %d\n", ndata, dim, kk);
    printf("+----------------------+\n\n");

    printf("\nGenerating random data...\n");
    for (int i = 0; i < ndata * dim; i++) {
        data[i] = (double)rand() / RAND_MAX;
        // printf("%f\t", data[i]);
        // if ((i + 1) % dim == 0) {
        //     printf("\n");
        // }
    }

    for (int i = 0; i < ndata; i++) {
        cluster_assign[i] = -1;
    }

    printf("\nCreating initial cluster centers...\n");
    initial_centers(dim, ndata, data, kk, &cluster_centroid);

    printf("\nSorting data with k-means...\n");
    error = kmeans(dim, ndata, &data, kk, &cluster_assign, &cluster_start, &cluster_size, &cluster_radius, &cluster_centroid);

    printf("\nSum of square errors: %f\n", error);

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
