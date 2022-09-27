#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int bipartition(int dim, int i0, int im, double* data, int chosen_dim, int cluster_id,      // input
    int cluster_start[2], int cluster_size[2],                                              // output
    double* cluster_bdry[2], double cluster_centroid,                                       // output
    short* cluster_assign) {                                                                // buffer

    // Buffer becomes the bi-partitioned data array for current two clusters.
    double* buffer = malloc(sizeof(double) * im * dim);

    printf("\ni0: %d\tim: %d\tdim: %d\tchosen_dim: %d\tcentroid: %f\tcluster ids: %d, %d\n\n", i0, im, dim, chosen_dim, cluster_centroid, cluster_id, (cluster_id + 1));

    int ii = 0,
        jj = im * dim,
        start_assign = i0 / dim,
        end_assign = (i0 / dim) + im - 1,
        is_left_first = dim,
        is_right_first = dim;

    for (int i = i0 + chosen_dim; i < im * dim + i0; i += dim) {
        if (data[i] <= cluster_centroid) {

            int left_bdry_idx = 0;
            for (int j = i - chosen_dim; j < i - chosen_dim + dim; j++) {
                // Update left cluster min and max boundaries in first occurrence.
                if (is_left_first > 0) {
                    cluster_bdry[cluster_id][left_bdry_idx] = data[j];
                    cluster_bdry[cluster_id][left_bdry_idx + 1] = data[j];
                    is_left_first -= 1;
                }
                // Updating left cluster min boundaries.
                else if (data[j] < cluster_bdry[cluster_id][left_bdry_idx]) {
                    cluster_bdry[cluster_id][left_bdry_idx] = data[j];
                }
                // Updating left cluster max boundaries.
                else if (data[j] > cluster_bdry[cluster_id][left_bdry_idx + 1]) {
                    cluster_bdry[cluster_id][left_bdry_idx + 1] = data[j];
                }

                buffer[ii] = data[j];
                left_bdry_idx += 2;
                ii += 1;
            }

            // Assigning current data point to left cluster.
            cluster_assign[start_assign] = cluster_id;
            start_assign += 1;
        }
        else {

            int right_bdry_idx = (dim * 2) - 1;
            for (int j = i - chosen_dim + dim - 1; j >= i - chosen_dim; j--) {
                // Update left cluster min and max boundaries in first occurrence.
                if (is_right_first > 0) {
                    cluster_bdry[cluster_id + 1][right_bdry_idx - 1] = data[j];
                    cluster_bdry[cluster_id + 1][right_bdry_idx] = data[j];
                    is_right_first -= 1;
                }
                // Updating right cluster min boundaries.
                else if (data[j] < cluster_bdry[cluster_id + 1][right_bdry_idx - 1]) {
                    cluster_bdry[cluster_id + 1][right_bdry_idx - 1] = data[j];
                }
                // Updating right cluster max boundaries.
                else if (data[j] > cluster_bdry[cluster_id + 1][right_bdry_idx]) {
                    cluster_bdry[cluster_id + 1][right_bdry_idx] = data[j];
                }

                buffer[jj - 1] = data[j];
                right_bdry_idx -= 2;
                jj -= 1;
            }

            // Assigning current data point to right cluster.
            cluster_assign[end_assign] = cluster_id + 1;
            end_assign -= 1;
        }
    }

    cluster_start[cluster_id] = i0;
    cluster_size[cluster_id] = ii / dim;

    cluster_start[cluster_id + 1] = i0 + ii;
    cluster_size[cluster_id + 1] = im - (ii / dim);

    ii = 0;
    // Updating original data array with cluster sorted data points
    for (int i = i0; i < im * dim + i0; i++) {
        data[i] = buffer[ii];
        ii += 1;
    }

    free(buffer);

    return 0;
}

int kd_tree(int dim, int ndata, double* data, int kk) {
    /***********************************************************************************
    dim -- number of dimensions or attributes of each datum.
    ndata -- total number of data points
    kk -- the number of subsets or clusters the whole dataset is partitioned into.
    Sizes of arrays:
        data[ndata*dim] -- the array storing ndata data points, each of dim dimensions.
        cluster_assign[ndata] -- array store index of the cluster a datum is assigned to.
        cluster_start[kk],
        cluster_size[kk],
        cluster_bdry[kk][2*dim],
        cluster_centroid[kk][dim]
    *************************************************************************************/

    int* cluster_size = malloc(sizeof(int) * kk);
    int* cluster_start = malloc(sizeof(int) * kk);

    short* cluster_assign = malloc(sizeof(short) * ndata);

    double** cluster_bdry = malloc(sizeof(double*) * kk * 2 * dim);
    double** cluster_centroid = malloc(sizeof(double*) * kk * dim);

    // Populating cluster_start & cluster_size
    cluster_start[0] = 0, cluster_size[0] = ndata;

    // Populating cluster_bdry
    for (int i = 0; i < kk; i++) {
        double* temp_bdry = malloc(sizeof(double) * dim * 2);
        for (int j = 0; j < dim * 2; j++) {
            if (j % 2 == 0) {
                temp_bdry[j] = __DBL_MAX__;
            }
            else {
                temp_bdry[j] = __DBL_MIN__;
            }
        }
        cluster_bdry[i] = temp_bdry;
    }

    double* centroid = malloc(sizeof(double) * dim);        // Temporarily stores centroid for current cluster.
    double* buffer = malloc(sizeof(double) * ndata);        // Temporarily stores data points from one dimension.

    int jj = 1;                                             // The number of clusters the dataset is partitioned into.

    while (jj < kk) {                                       // At the beginning, jj=1, meaning only one cluster, the whole dataset.

        int* c_start_buffer = malloc(sizeof(int) * kk);
        int* c_size_buffer = malloc(sizeof(int) * kk);

        for (int j = 0; j < jj; j++) {                      // j loops through indices of all jj clusters.

            int chosen_dim = -1;
            double largest_variance = __DBL_MIN__;

            for (int i = 0; i < dim; i++) {

                double sum = 0.0;
                int k = 0;

                for (int ii = cluster_start[j] + i; ii < cluster_start[j] + cluster_size[j] * dim; ii += dim) {
                    buffer[k] = data[ii];
                    sum += data[ii];
                    k += 1;
                }

                // Mean (centroid) calculation for current dimension
                centroid[i] = sum / k;

                double variance_numerator = 0.0;
                for (int ii = 0; ii < cluster_size[j]; ii++) {
                    variance_numerator += pow((buffer[ii] - centroid[i]), 2);
                }

                // Variance calculation for current dimension
                double dim_variance = variance_numerator / k;

                printf("Dimension %d --> Variance: %f, Mean: %f\n", i, dim_variance, centroid[i]);

                // Determines dimension with largest variance
                if (dim_variance > largest_variance) {
                    largest_variance = dim_variance;
                    chosen_dim = i;
                }
            }

            cluster_centroid[j] = centroid;

            int cluster_id = (j == 0) ? j : j + 1;

            //   call bipartition() to partition the j-th cluster into 2 clusters
            bipartition(dim, cluster_start[j], cluster_size[j], data, chosen_dim, cluster_id,
                c_start_buffer, c_size_buffer, cluster_bdry, cluster_centroid[j][chosen_dim], cluster_assign);


            printf("Bipartition result:\n");
            for (int i = 0; i < ndata * dim; i++) {
                printf("%f\t", data[i]);
                if ((i + 1) % dim == 0) {
                    printf("\n");
                }
            }
            printf("\n");

        }

        cluster_start = c_start_buffer;
        cluster_size = c_size_buffer;

        printf("\n");

        for (int i = 0; i < kk; i++) {
            printf("Cluster %d --> Start: %d\tSize: %d\n", i, cluster_start[i], cluster_size[i]);
        }

        printf("\n");

        // for (int i = 0; i < jj * 2; i++) {
        //     for (int j = 0; j < dim * 2; j += 2) {
        //         printf("Cluster %d / Dimension %d --> Min %f, Max: %f\n", i, j / 2, cluster_bdry[i][j], cluster_bdry[i][j + 1]);
        //     }
        // }

        for (int i = 0; i < ndata; i++) {
            printf("%d\n", cluster_assign[i]);
        }

        printf("\n+-----------------------------------------------+\n\n");

        jj *= 2;
    }

    printf("\n");

    free(buffer);
    free(centroid);
    free(cluster_bdry);
    free(cluster_size);
    free(cluster_start);
    free(cluster_centroid);

    return 0;
}

// int search_kdtree(int dim, int ndata, double *data, int kk,
//                   int *cluster_start, int *cluster_size, double **cluster_bdry,
//                   double *query_pt, double *result_pt) {
// 
//     /* search_kdtree() returns the number of data points checked */
// }

int main() {

    int ndata = 10, dim = 3, kk = 4;
    double* data = malloc(sizeof(double) * ndata * dim);

    // Generating data
    printf("\n+------------------------- Original data -------------------------+\n");
    for (int i = 0; i < ndata * dim; i++) {
        data[i] = (double)rand() / RAND_MAX;
        printf("%f\t", data[i]);

        if ((i + 1) % dim == 0) {
            printf("\n");
        }
    }
    printf("+-----------------------------------------------------------------+\n\n");


    kd_tree(dim, ndata, data, kk);

    free(data);

    return 0;
}
