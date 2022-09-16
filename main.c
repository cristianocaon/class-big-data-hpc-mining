#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int bipartition(int dim, int i0, int im, double *data, int chosen_dim, //  input
                int cluster_start[2], int cluster_size[2],             // output
                double *cluster_bdry[2], double *cluster_centroid[2],  // output
                short *cluster_assign)                                 // buffer
{

    return 0;
}

int kd_tree(int dim, int ndata, double *data, int kk,
            int *cluster_start, int *cluster_size, double **cluster_bdry, double **cluster_centroid,
            short *cluster_assign)
{
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

    double *centroids = malloc(sizeof(double) * dim);
    double *buffer = malloc(sizeof(double) * ndata); // Temporarily stores data points from one dimension

    int jj = 1; // the number of clusters the dataset is partitioned into.

    while (jj < kk) // At the beginning, jj=1, meaning only one cluster, the whole dataset.
    {
        for (int j = 0; j < jj; j++) // j loops through indices of all jj clusters
        {
            int chosen_dim = -1;
            double largest_variance = -1.0, chosen_centroid = -1.0;

            for (int i = 0; i < dim; i++)
            {
                double sum = 0.0;
                int k = 0;
                for (int ii = cluster_start[j] + i; ii < cluster_size[j] * dim; ii += dim)
                {
                    buffer[k] = data[ii];
                    sum += data[ii];
                    k += 1;
                }

                // Mean (centroid) calculation for current dimension
                centroids[i] = sum / k;

                double variance_numerator = 0.0;
                for (int ii = 0; ii < cluster_size[j]; ii++)
                {
                    variance_numerator += pow((buffer[ii] - centroids[i]), 2);
                }

                // Variance calculation for current dimension
                double dim_variance = variance_numerator / k;

                printf("Current dimension: %d\tVariance: %f\tMean: %f\n", i, dim_variance, centroids[i]);

                // Determines dimension with largest variance
                if (dim_variance > largest_variance)
                {
                    largest_variance = dim_variance;
                    chosen_centroid = centroids[i];
                    chosen_dim = i;
                }
            }

            printf("Largest variance is from dimension %d: %f\n\n", chosen_dim, largest_variance);
            //   call bipartition() to partition the j-th cluster into 2 clusters
        }
        jj *= 2;
    }

    return 0;
}

// int search_kdtree(int dim, int ndata, double *data, int kk,
//                   int *cluster_start, int *cluster_size, double **cluster_bdry,
//                   double *query_pt, double *result_pt)
// {
//     /* search_kdtree() returns the number of data points checked */
// }

int main()
{
    int ndata = 1000000, dim = 16, kk = 3;

    double *data = malloc(sizeof(double) * ndata * dim);

    int *cluster_start = malloc(sizeof(int) * kk);
    int *cluster_size = malloc(sizeof(int) * kk);

    double **cluster_bdry = malloc(sizeof(double *) * kk * 2 * dim);
    double **cluster_centroid = malloc(sizeof(double *) * kk * dim);

    short *cluster_assign = malloc(sizeof(short) * ndata);

    // Generating data
    for (int i = 0; i < ndata * dim; i++)
    {
        data[i] = (double)rand() / RAND_MAX;
    }

    // Populating cluster_start & cluster_size
    cluster_start[0] = 0, cluster_size[0] = ndata;
    for (int i = 1; i < kk; i++)
    {
        cluster_start[i] = -1;
        cluster_size[i] = -1;
    }

    // Populating cluster_assign
    for (int i = 0; i < ndata; i++)
    {
        cluster_assign[i] = 0;
    }

    kd_tree(dim, ndata, data, kk, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_assign);

    return 0;
}
