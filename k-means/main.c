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


int initial_centers(int dim, int i0, int im, double* data, int k,
    double** cluster_centroid) {
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

    return 0;
}


double kmeans(int dim, int i0, int im, double* data, int k, // input
    short* cluster_assign,                                  // buffer
    int* cluster_start, int* cluster_size,                  // output
    double* cluster_radius, double** cluster_centroid) {    // output
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

    return 0;
}


int search_kmeans(int dim, int ndata, double* data, int k,
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

    return 0;
}
