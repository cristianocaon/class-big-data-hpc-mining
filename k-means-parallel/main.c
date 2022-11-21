/******************************************

 * Texas Tech University
 * CS 5331-001: Mining Big Data with HPC
 * Professor: Yu Zhuang
 * Assignment 4: Parallel K-Means
 * Student: Cristiano Eleutherio Caon
 * R#: 11474435

******************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// Data structure to send a value and process rank
struct double_rank {
  double value;
  int    rank;
};


int initial_centers(int dim, int ndata, double* data, int kk,
  double*** cluster_centroid, int rank, int np) {
  /******************************************************
   * Returns the initial centers for the 'kk' clusters.
  ******************************************************/

  int root = 0, largest_idx = 0, clusters = 2;

  double temp_dist;
  double* d_min = (double*)malloc(sizeof(double) * ndata);

  struct double_rank local_dist;
  struct double_rank global_dist;

  // Generating first cluster centroid
  (*cluster_centroid)[0] = (double*)malloc(sizeof(double) * dim);
  if (rank == root) {
    for (int i = 0; i < dim; i++) {
      (*cluster_centroid)[0][i] = data[i];
    }
  }

  // Broadcasting first cluster center
  MPI_Bcast((*cluster_centroid)[0], dim, MPI_DOUBLE, root, MPI_COMM_WORLD);

  local_dist.value = __DBL_MIN__;
  local_dist.rank = rank;
  // Finding second furthest away centroid
  (*cluster_centroid)[1] = (double*)malloc(sizeof(double) * dim);
  for (int i = 0; i < ndata * dim; i += dim) {
    temp_dist = 0.0;
    for (int j = i, k = 0; j < i + dim; j++, k++) {
      temp_dist += pow(fabs((*cluster_centroid)[0][k]) - fabs(data[j]), 2);
    }
    temp_dist = sqrt(temp_dist);
    if (temp_dist > local_dist.value) {
      for (int j = i, k = 0; j < i + dim; j++, k++) {
        (*cluster_centroid)[1][k] = data[j];
      }
      local_dist.value = temp_dist;
    }
  }

  // Obtaining global MAX distance
  MPI_Allreduce(&local_dist, &global_dist, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  // Broadcasting second furthest away centroid
  MPI_Bcast((*cluster_centroid)[1], dim, MPI_DOUBLE, global_dist.rank, MPI_COMM_WORLD);

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
    // Getting furthest data point to a center distance
    largest_idx = 0;
    local_dist.value = d_min[0];
    for (int i = 1; i < ndata; i++) {
      if (d_min[i] > local_dist.value) {
        local_dist.value = d_min[i];
        largest_idx = i;
      }
    }

    // Obtaining nth global MAX distance
    MPI_Allreduce(&local_dist, &global_dist, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    // Saving data point as a new cluster center
    (*cluster_centroid)[clusters] = (double*)malloc(sizeof(double) * dim);
    if (rank == global_dist.rank) {
      for (int i = largest_idx * dim, j = 0; i < largest_idx * dim + dim; i++, j++) {
        (*cluster_centroid)[clusters][j] = data[i];
      }
    }

    // Broadcasting nth furthest away centroid
    MPI_Bcast((*cluster_centroid)[clusters], dim, MPI_DOUBLE, global_dist.rank, MPI_COMM_WORLD);

    clusters++;
  }

  free(d_min);

  return 0;
}


double kmeans(int dim, int ndata, double** data, int kk, // input
  short** cluster_assign,                                // buffer
  int** cluster_start, int** cluster_size,               // output
  double** cluster_radius, double*** cluster_centroid,   // output
  int rank, int np) {
  /****************************************************
   * Organizes the data array using K-Means algorithm
   *  into 'kk' number of clusters.
   *
   * Returns the sum of square errors.
  ****************************************************/

  int all_proc_cluster_change, all_proc_count;
  int count, count_cluster_change, chosen_cluster, stop_iteration = 0;

  double all_proc_total_error;
  double temp_dist, d_min, total_error, cluster_error, point_error;

  double* data_buffer;
  double* centroid_buffer;
  double* all_proc_centroid_buffer;

  struct double_rank local_radius;
  struct double_rank global_radius;

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

    // Getting max count cluster change from processes
    MPI_Allreduce(&count_cluster_change, &all_proc_cluster_change, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // Exit condition
    if (all_proc_cluster_change == 0) {
      stop_iteration = 1;
      break;
    }

    // Re-calculating cluster centers
    for (int i = 0; i < kk; i++) {
      all_proc_count = 0;
      count = 0;
      all_proc_centroid_buffer = (double*)malloc(sizeof(double) * dim);
      centroid_buffer = (double*)malloc(sizeof(double) * dim);
      for (int j = 0, k = 0; j < ndata * dim; j += dim, k++) {
        if ((*cluster_assign)[k] == i) {
          for (int ii = 0; ii < dim; ii++) {
            centroid_buffer[ii] += (*data)[j + ii];
          }
          count++;
        }
      }

      // Adding all dimensions from all processes
      MPI_Allreduce(centroid_buffer, all_proc_centroid_buffer, dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      // Adding count change from all processes
      MPI_Allreduce(&count, &all_proc_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      // Updating cluster centroid with info from all processes
      for (int j = 0; j < dim; j++) {
        all_proc_centroid_buffer[j] /= all_proc_count;
        (*cluster_centroid)[i][j] = all_proc_centroid_buffer[j];
      }
    }

    // Re-calculating cluster radius
    for (int i = 0; i < kk; i++) {
      local_radius.value = __DBL_MIN__;
      local_radius.rank = rank;
      for (int j = 0, k = 0; j < ndata * dim; j += dim, k++) {
        if ((*cluster_assign)[k] == i) {
          temp_dist = 0.0;
          for (int ii = 0; ii < dim; ii++) {
            temp_dist += pow(fabs((*cluster_centroid)[i][ii]) - fabs((*data)[j + ii]), 2);
          }
          temp_dist = sqrt(temp_dist);
          local_radius.value = (temp_dist > local_radius.value) ? temp_dist : local_radius.value;
        }
      }
      // Getting largest radius from processes
      MPI_Allreduce(&local_radius, &global_radius, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      // Broadcasting largest radius to all processes
      MPI_Bcast(&local_radius.value, 1, MPI_DOUBLE, local_radius.rank, MPI_COMM_WORLD);
      (*cluster_radius)[i] = local_radius.value;
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
        point_error += fabs((*cluster_centroid)[i][ii]) - fabs((*data)[k]);
      }
      if (point_error > 0) {
        cluster_error += pow(point_error, 2);
      }
    }
    if (cluster_error > 0) {
      total_error += cluster_error / (*cluster_size)[i];
    }
  }

  // Getting all processes total error
  MPI_Allreduce(&total_error, &all_proc_total_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  free(centroid_buffer);
  free(all_proc_centroid_buffer);

  return all_proc_total_error;
}


int search_kmeans(int dim, int ndata, double* data, int kk,
  int* cluster_start, int* cluster_size,
  double* cluster_radius, double** cluster_centroid,
  double* query_pt, double** result_pt,
  int rank, int np) {
  /*****************************************************
   * Searches for closest data point from data array to
   *  given query points.
   *
   * Returns the number of data points checked
  *****************************************************/

  int global_checked = 0, local_checked = 0, cluster = -1;

  double radius, temp_dist;

  struct double_rank local_dist;
  struct double_rank global_dist;

  local_dist.rank = rank;
  local_dist.value = __DBL_MAX__;

  // Finding closest cluster center
  for (int i = 0; i < kk; i++) {
    temp_dist = 0.0;
    for (int j = 0; j < dim; j++) {
      temp_dist += pow(fabs(query_pt[j]) - fabs(cluster_centroid[i][j]), 2);
    }
    temp_dist = sqrt(temp_dist);
    if (temp_dist < local_dist.value) {
      local_dist.value = temp_dist;
      cluster = i;
    }
  }

  // Finding closest data point inside cluster
  local_dist.value = __DBL_MAX__;
  for (int i = cluster_start[cluster]; i < cluster_start[cluster] + cluster_size[cluster] * dim; i += dim) {
    temp_dist = 0.0;
    for (int j = i, k = 0; j < i + dim; j++, k++) {
      temp_dist += pow(fabs(query_pt[k]) - fabs(data[j]), 2);
    }
    temp_dist = sqrt(temp_dist);
    if (temp_dist < local_dist.value) {
      local_dist.value = temp_dist;
      for (int j = i, k = 0; j < i + dim; j++, k++) {
        (*result_pt)[k] = data[j];
      }
    }
    local_checked++;
  }

  radius = local_dist.value;

  // Looking for even closer data point within query point cluster
  for (int i = 0; i < kk; i++) {
    temp_dist = 0.0;
    for (int j = 0; j < dim; j++) {
      temp_dist += pow(fabs(query_pt[j]) - fabs(cluster_centroid[i][j]), 2);
    }
    temp_dist = sqrt(temp_dist);
    // Check if query point cluster intersects with other clusters
    if (temp_dist < radius + cluster_radius[i]) {
      for (int j = cluster_start[i]; j < cluster_start[i] + cluster_size[i] * dim; j += dim) {
        temp_dist = 0.0;
        for (int k = j, ii = 0; k < j + dim; k++, ii++) {
          temp_dist += pow(fabs(query_pt[ii]) - fabs(data[k]), 2);
        }
        temp_dist = sqrt(temp_dist);
        // Checking if there is a closer data point
        if (temp_dist < local_dist.value) {
          local_dist.value = temp_dist;
          for (int k = j, ii = 0; k < j + dim; k++, ii++) {
            (*result_pt)[ii] = data[k];
          }
        }
        local_checked++;
      }
    }
  }

  // Getting global closest distance
  MPI_Allreduce(&local_dist, &global_dist, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

  // Broadcasting result point of closest distance
  if (rank == global_dist.rank) {
    for (int i = 0; i < np; i++) {
      if (i != global_dist.rank) {
        MPI_Send((*result_pt), dim, MPI_DOUBLE, i, 900, MPI_COMM_WORLD);
      }
    }
  }
  else {
    // Receiving closest point
    MPI_Recv((*result_pt), dim, MPI_DOUBLE, global_dist.rank, 900, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  // Getting total number of data points checked from all processes
  MPI_Allreduce(&local_checked, &global_checked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  return global_checked;
}


int main(int argc, char** argv) {
  /***********************************************
   * Steps performed in main program:
   *
   *  1- Generate random data
   *  2- Create initial cluster centers
   *  3- Sort data with K-Means
   *  4- Search closest distance to query points
  ***********************************************/

  int np, rank;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  int seed = 0, checked = 0, ndata = 100000, dim = 4, kk = 10;

  int* cluster_size = (int*)malloc(sizeof(int) * kk);
  int* cluster_start = (int*)malloc(sizeof(int) * kk);

  double error;

  double* data = (double*)malloc(sizeof(double) * (ndata / np) * dim);
  double* query_pt = (double*)malloc(sizeof(double) * dim);
  double* result_pt = (double*)malloc(sizeof(double) * dim);
  double* cluster_radius = (double*)malloc(sizeof(double) * kk);
  double** cluster_centroid = (double**)malloc(sizeof(double*) * kk);

  short* cluster_assign = (short*)malloc(sizeof(short) * ndata / np);

  if (rank == 0) {
    printf("\nParameters:\n+----------------------+");
    printf("\nndata = %d\ndim = %d\nkk = %d\n", ndata, dim, kk);
    printf("+----------------------+\n");
    printf("\nGenerating random data...");
  }

  srand(rank + 1);
  for (int i = 0; i < (ndata / np) * dim; i++) {
    data[i] = (double)rand() / RAND_MAX;
  }

  if (rank == 0) {
    printf("\n\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n");
    printf("\nCreating initial cluster centers...");
  }
  initial_centers(dim, ndata / np, data, kk, &cluster_centroid, rank, np);

  for (int i = 0; i < ndata / np; i++) {
    cluster_assign[i] = -1;
  }

  if (rank == 0) {
    printf("\nSorting data with k-means...");
  }
  error = kmeans(dim, ndata / np, &data, kk, &cluster_assign, &cluster_start,
    &cluster_size, &cluster_radius, &cluster_centroid, rank, np);

  if (rank == 0) {
    printf("\nSum of square errors: %f", error);
    printf("\n\n+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n");
    printf("\nGenerating a random query point...");

    // Broadcasting seed for query point generation
    seed = np * np;
    for (int i = 1; i < np; i++) {
      MPI_Send(&seed, 1, MPI_INT, i, 999, MPI_COMM_WORLD);
    }
  }
  else {
    // Receiving seed from root process for query point
    MPI_Recv(&seed, 1, MPI_INT, 0, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  srand(seed);
  for (int i = 0; i < dim; i++) {
    query_pt[i] = (double)rand() / RAND_MAX;
  }

  if (rank == 0) {
    printf("\nSearching for closest data point to query point...");
  }
  checked = search_kmeans(dim, ndata / np, data, kk, cluster_start, cluster_size,
    cluster_radius, cluster_centroid, query_pt, &result_pt, rank, np);

  if (rank == 0) {
    printf("\nChecked %d data points for query point.\n", checked);
    printf("\nResult point: \n");
    for (int i = 0; i < dim; i++) {
      printf("Coordinate %d\t --> %f\n", i, result_pt[i]);
    }
  }

  // Clean up
  for (int i = 0; i < kk; i++) {
    free(cluster_centroid[i]);
  }
  free(data);
  free(query_pt);
  free(result_pt);
  free(cluster_size);
  free(cluster_start);
  free(cluster_assign);
  free(cluster_radius);

  MPI_Finalize();

  return 0;
}
