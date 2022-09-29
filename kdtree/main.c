/******************************************

 * Texas Tech University
 * CS 5331-001: Mining Big Data with HPC
 * Professor: Yu Zhuang
 * Assignment 1: K-D tree
 * Student: Cristiano Eleutherio Caon
 * R#: 11474435

******************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int bipartition(int dim, int i0, int im, double* data, int chosen_dim, int cluster_id, int cluster_start[2],
    int cluster_size[2], double* cluster_bdry[2], double cluster_centroid, short* cluster_assign) {
    /**********************************************
     * Partitions given cluster into two clusters.
    **********************************************/

    double* partitioned_data = malloc(sizeof(double) * im * dim);

    int left_bdry_idx,
        right_bdry_idx,
        ii = 0,
        jj = im * dim,
        start_assign = i0 / dim,
        end_assign = (i0 / dim) + im - 1,
        is_left_first = dim,
        is_right_first = dim;

    // Loop through cluster data points on chosen dimension.
    for (int i = i0 + chosen_dim; i < im * dim + i0; i += dim) {
        // Check assignment to left-cluster.
        if (data[i] <= cluster_centroid) {
            left_bdry_idx = 0;
            // Loop forwards through current data point.
            for (int j = i - chosen_dim; j < i - chosen_dim + dim; j++) {
                // Update left-cluster min and max boundaries in first occurrence.
                if (is_left_first > 0) {
                    cluster_bdry[cluster_id][left_bdry_idx] = data[j];
                    cluster_bdry[cluster_id][left_bdry_idx + 1] = data[j];
                    is_left_first -= 1;
                }
                // Check left-cluster min boundaries.
                else if (data[j] < cluster_bdry[cluster_id][left_bdry_idx]) {
                    cluster_bdry[cluster_id][left_bdry_idx] = data[j];
                }
                // Check left-cluster max boundaries.
                else if (data[j] > cluster_bdry[cluster_id][left_bdry_idx + 1]) {
                    cluster_bdry[cluster_id][left_bdry_idx + 1] = data[j];
                }
                partitioned_data[ii] = data[j];
                left_bdry_idx += 2;
                ii += 1;
            }
            // Assign current data point to left-cluster.
            cluster_assign[start_assign] = cluster_id;
            start_assign += 1;
        }
        // Assignment to right-cluster.
        else {
            right_bdry_idx = (dim * 2) - 1;
            // Loop backwards through current data point.
            for (int j = i - chosen_dim + dim - 1; j >= i - chosen_dim; j--) {
                // Update right-cluster min and max boundaries in first occurrence.
                if (is_right_first > 0) {
                    cluster_bdry[cluster_id + 1][right_bdry_idx - 1] = data[j];
                    cluster_bdry[cluster_id + 1][right_bdry_idx] = data[j];
                    is_right_first -= 1;
                }
                // Check right-cluster min boundaries.
                else if (data[j] < cluster_bdry[cluster_id + 1][right_bdry_idx - 1]) {
                    cluster_bdry[cluster_id + 1][right_bdry_idx - 1] = data[j];
                }
                // Check right-cluster max boundaries.
                else if (data[j] > cluster_bdry[cluster_id + 1][right_bdry_idx]) {
                    cluster_bdry[cluster_id + 1][right_bdry_idx] = data[j];
                }
                partitioned_data[jj - 1] = data[j];
                right_bdry_idx -= 2;
                jj -= 1;
            }
            // Assign current data point to right-cluster.
            cluster_assign[end_assign] = cluster_id + 1;
            end_assign -= 1;
        }
    }

    // Update left-cluster start and size.
    cluster_start[cluster_id] = i0;
    cluster_size[cluster_id] = ii / dim;
    // Update right-cluster start and size.
    cluster_start[cluster_id + 1] = i0 + ii;
    cluster_size[cluster_id + 1] = im - (ii / dim);

    ii = 0;

    // Update original data array.
    for (int i = i0; i < im * dim + i0; i++) {
        data[i] = partitioned_data[ii];
        ii += 1;
    }

    free(partitioned_data);

    return 0;
}

int kd_tree(int dim, int ndata, double* data, int kk, int* cluster_start, int* cluster_size,
    double** cluster_bdry, double** cluster_centroid, short* cluster_assign) {
    /********************************************
     * Sorts a given multidimensional data array
       using K-D tree.
    ********************************************/

    int cluster_id, chosen_dim, k, jj = 1;

    double largest_variance, sum, variance_numerator, dim_variance;

    int* size_buffer = malloc(sizeof(int) * kk);
    int* start_buffer = malloc(sizeof(int) * kk);

    double* temp_bdry;
    double* centroid = malloc(sizeof(double) * dim);
    double* dim_pts = malloc(sizeof(double) * ndata);

    // Populating first cluster start and size.
    cluster_start[0] = 0, cluster_size[0] = ndata;

    // Populating clusters boundaries.
    for (int i = 0; i < kk; i++) {
        temp_bdry = malloc(sizeof(double) * dim * 2);
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

    // Loop until desired number of clusters.
    while (jj < kk) {
        cluster_id = 0;
        // Loop through partitioned clusters.
        for (int j = 0; j < jj; j++) {
            chosen_dim = -1;
            largest_variance = __DBL_MIN__;
            // Loop through dimensions.
            for (int i = 0; i < dim; i++) {
                k = 0;
                sum = 0.0;
                // Loop through current cluster dimension coordinates to calculate mean.
                for (int ii = cluster_start[j] + i; ii < cluster_start[j] + cluster_size[j] * dim; ii += dim) {
                    dim_pts[k] = data[ii];
                    sum += data[ii];
                    k += 1;
                }

                // Calculates the Mean (centroid) for current cluster dimension.
                centroid[i] = sum / k;
                variance_numerator = 0.0;

                // Loop through current cluster dimension coordinates to calculate variance.
                for (int ii = 0; ii < cluster_size[j]; ii++) {
                    variance_numerator += pow((dim_pts[ii] - centroid[i]), 2);
                }

                // Variance calculation for current cluster dimension.
                dim_variance = variance_numerator / k;

                // Check for dimension with largest variance.
                if (dim_variance > largest_variance) {
                    largest_variance = dim_variance;
                    chosen_dim = i;
                }
            }

            cluster_centroid[j] = centroid;

            // Bipartition the j-th cluster into 2 clusters
            bipartition(dim, cluster_start[j], cluster_size[j], data, chosen_dim,
                cluster_id, start_buffer, size_buffer, cluster_bdry,
                cluster_centroid[j][chosen_dim], cluster_assign);

            cluster_id += 2;
        }

        // Update clusters start and size.
        for (int i = 0; i < kk; i++) {
            cluster_start[i] = start_buffer[i];
            cluster_size[i] = size_buffer[i];
        }

        jj *= 2;
    }

    free(dim_pts);
    free(centroid);
    free(size_buffer);
    free(start_buffer);

    return 0;
}


int search_kd_tree(int dim, int ndata, double* data, int kk, int* cluster_start, int* cluster_size,
    double** cluster_bdry, double* query_pts, double* result_pts, int num_query_pts) {
    /********************************************
     * Searches K-D tree for closest data points
       to given query points.
    ********************************************/

    int k, count_cluster_dims, distance, dist_dims;

    double pts_checked = 0.0;

    int* dim_bdry_flags;

    double* d_mins = malloc(sizeof(double) * num_query_pts);

    // Populating min distances for each query point.
    for (int ii = 0; ii < num_query_pts; ii++) {
        d_mins[ii] = __DBL_MAX__;
    }

    // Loop through query points.
    for (int ii = 0; ii < num_query_pts * dim; ii += dim) {
        // Loop through clusters.
        for (int i = 0; i < kk; i++) {
            k = 0;
            dist_dims = 0.0;
            count_cluster_dims = 0;
            dim_bdry_flags = malloc(sizeof(int) * dim);
            // Loop through dimensions boundaries.
            for (int j = 0; j < dim * 2; j += 2) {
                // Check if query dimension is within current cluster dimension boundaries.
                if ((query_pts[ii + k] >= cluster_bdry[i][j]) &&
                    (query_pts[ii + k] <= cluster_bdry[i][j + 1])) {
                    count_cluster_dims += 1;
                    dim_bdry_flags[k] = 0;
                }
                // If query dimension is less than min cluster boundary, get query dimension distance to min cluster boundary.
                else if (query_pts[ii + k] < cluster_bdry[i][j]) {
                    dim_bdry_flags[k] = 1;
                }
                // If query dimension is greater than max cluster boundary, get query dimension distance to max cluster boundary.
                else if (query_pts[ii + k] > cluster_bdry[i][j + 1]) {
                    dim_bdry_flags[k] = 2;
                }
                k += 1;
            }

            k = 0;

            // Check if query point is within all cluster dimensions boundaries, i.e., inside cluster.
            if (count_cluster_dims == dim) {
                // Loop through cluster data points.
                for (int j = cluster_start[i]; j < cluster_start[i] + cluster_size[i] * dim; j++) {

                    // Cluster data point & query point current dimension distance.
                    dist_dims += pow(fabs(data[j]) - fabs(query_pts[ii + k]), 2);

                    if (k == (dim - 1)) {
                        pts_checked += 1;
                        distance = sqrt(dist_dims);
                        // Check if data point distance is the new smallest.
                        if (distance < d_mins[ii / dim]) {
                            k = 0;
                            // Loop through data point of smallest distance to save coordinates.
                            for (int jj = j; jj < j + dim; jj++) {
                                result_pts[ii + k] = data[jj];
                                k += 1;
                            }
                            d_mins[ii / dim] = distance;
                            k = 0;
                        }
                    }
                    else {
                        k += 1;
                    }
                }
            }
            // Query point is outside of cluster boundaries.
            else {
                // Loop through dimensions
                for (int j = 0; j < dim * 2; j += 2) {
                    // Get query dimension distance to min cluster boundary.
                    if (dim_bdry_flags[k] == 1) {
                        dist_dims += pow(fabs(cluster_bdry[i][j]) - fabs(query_pts[ii + k]), 2);
                    }
                    // Get query dimension distance to max cluster boundary.
                    else if (dim_bdry_flags[k] == 2) {
                        dist_dims += pow(fabs(cluster_bdry[i][j + 1]) - fabs(query_pts[ii + k]), 2);
                    }
                    k += 1;
                }
                pts_checked += 1;
                distance = sqrt(dist_dims);

                // Check if distance to cluster is the new smallest.
                if (distance < d_mins[ii / dim]) {
                    d_mins[ii / dim] = distance;
                }
            }
        }
    }

    printf("\n\nResults:\n+----------------------------------------------------------------------------+");
    printf("\nTotal number of data points checked for %d query points: %.0f\n", num_query_pts, pts_checked);
    printf("\nAverage number of data points checked for %d query points: %.2f\n", num_query_pts, (pts_checked / num_query_pts));
    printf("\nPercentage of total data points checked: %.0f / (%d * %d) = %.2f%%\n", pts_checked, ndata, num_query_pts, pts_checked / (ndata * num_query_pts) * 100);
    printf("+----------------------------------------------------------------------------+\n");

    free(d_mins);
    free(dim_bdry_flags);

    return 0;
}


int main() {

    int ndata = 1000000, dim = 16, kk = 1024, num_query_pts = 1000;

    int* cluster_size = malloc(sizeof(int) * kk);
    int* cluster_start = malloc(sizeof(int) * kk);

    short* cluster_assign = malloc(sizeof(short) * ndata);

    double* data = malloc(sizeof(double) * ndata * dim);
    double* query_pts = malloc(sizeof(double) * dim * num_query_pts);
    double* result_pts = malloc(sizeof(double) * dim * num_query_pts);
    double** cluster_bdry = malloc(sizeof(double*) * kk * 2 * dim);
    double** cluster_centroid = malloc(sizeof(double*) * kk * dim);

    printf("\nParameters:\n+----------------------+");
    printf("\nndata = %d\ndim = %d\nkk = %d\nnum_query_pts = %d\n", ndata, dim, kk, num_query_pts);
    printf("+----------------------+\n\n");

    printf("\nGenerating random data...\n");
    for (int i = 0; i < ndata * dim; i++) {
        data[i] = (double)rand() / RAND_MAX;
    }

    printf("\nSorting data using K-D tree...\n");
    kd_tree(dim, ndata, data, kk, cluster_start, cluster_size, cluster_bdry,
        cluster_centroid, cluster_assign);

    printf("\nGenerating random query points...\n");
    for (int i = 0; i < num_query_pts * dim; i++) {
        query_pts[i] = (double)rand() / RAND_MAX;
    }

    printf("\nSearching for closest data points...\n");
    search_kd_tree(dim, ndata, data, kk, cluster_start, cluster_size,
        cluster_bdry, query_pts, result_pts, num_query_pts);

    free(data);
    free(query_pts);
    free(result_pts);
    free(cluster_bdry);
    free(cluster_size);
    free(cluster_start);
    free(cluster_assign);
    free(cluster_centroid);

    return 0;
}
