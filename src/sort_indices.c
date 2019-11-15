/**
 * Provide sorting functions (currently CPU). 
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sort_indices.h>

// dim = 2, visibility limit to compilation unit
static const int dim = 2;

void swapf(float* a, float* b) {
  float t = *a;
  *a = *b;
  *b = t;
}

void swapi(int* a, int* b) {
  int t = *a;
  *a = *b;
  *b = t;
}

/* This function takes last element as pivot, places the pivot element at its correct position in sorted
 * array, and places all smaller (smaller than pivot) to left of pivot and all greater elements to right
 * of pivot.
 */
int partition(float *values, int *indices, int low, int high)
{
  float pivot = values[high];
  // i = index of smaller element
  int i = (low - 1);  

  for (int j = low; j <= high- 1; j++)
  {
    if (values[j] <= pivot)
    {
      i++;
      swapf(&values[i], &values[j]);
      swapi(&indices[i], &indices[j]);
    }
  }
  i++;
  swapf(&values[i], &values[high]);
  swapi(&indices[i], &indices[high]);
  return i;
}

/* The main function that implements quick sort. 
*/
void quicksort(float *values, int *indices, int low, int high)
{
  if (low < high)
  {
    // pi is partitioning index, values[p] is now at right place 
    int pi = partition(values, indices, low, high);

    // Separately sort elements before partition and after partition
    quicksort(values, indices, low, pi - 1);
    quicksort(values, indices, pi + 1, high);
  }
}

/* Tensorflow has also argsort with slightly different arguments.
 *
 * This uses memcpy because sorting is in-place. Hopefully tf.argsort is implemented differently.
 */
void argsort(int n, const float* values, int *indices) {
  // set indices to a consecutive sequence from 0 to n-1
  for (int i = 0; i < n; ++i) {
    indices[i] = i;
  }
  float *tmp = (float*)malloc(sizeof(float)*n);
  for (int i = 0; i < n; ++i) {
    tmp[i] = values[i];
  }
  quicksort(tmp, indices, 0, n-1);
  free(tmp);
}

/**
 * Number of values smaller than threshold.
 * Assumes increasing order of values[indices[.]].
 * Note, we return the index itself, not the value within the array of indices. 
 */
int count_items_below_threshold(int n, const float* values, const int *indices, float threshold) {
  int result = 0;
  for (int i = 0; i < n; ++i) {
    if (values[indices[i]] > threshold) {
      break;
    } 
    result++;
  }
  return result;
}

/**
 * Average. We can ignore indices here.
 */
float avg(int n, const float* values, int dim_index) {
  long double mean = 0;
  for (int i = 0; i < n; ++i) {
    mean += values[i * dim + dim_index];
  }
  float result = mean/n;
  printf("Average [dim=%i]: %f\n", dim_index, result);
  return result;
}

/**
 * Smoothing average where we calculate a window of size k, with k nearest neighbours. 
 */
void smoothing_avg(int n, const float* values, const int *indices, float *mean, int k) {
  if (k == 0) {
    for (int i = 0; i < n; ++i) {
      mean[i] = values[i];
    }
    return;
  }
//  assert (k < n);
//  assert (k % 2 == 0);

  for (int i = 0; i < n; ++i) {
    mean[i] = 0;
  }

  int k2 = k/2;
  // first sum for first entry up to k/2 + 1 values
  for (int i = 0; i < k2 + 1; ++i) {
    mean[indices[0]] += values[indices[i]];
  }
  // then from second entry on calculate moving average by adding item k/2 forwards
  // normalize last item considered i-1 by size of (still increasing sliding window)
  for (int i = 1; i < k2 + 1; ++i) {
    mean[indices[i]] = mean[indices[i-1]] + values[indices[i + k2]];
    mean[indices[i-1]] /= (k2 + i);
  }
  // from k/2 + 1 on our sliding window is full size (k + 1), we subtract item -(k/2+1) and add item k/2
  // normalize last item with this window size
  for (int i = k2 + 1; i < n - k2; ++i) {
    mean[indices[i]] = mean[indices[i-1]] + values[indices[i + k2]] - values[indices[i - k2 - 1]];
    mean[indices[i-1]] /= (k + 1);
  }
  // now are sliding window is decreasing in size again, only subtract items 
  for (int i = n - k2; i < n; ++i) {
    mean[indices[i]] = mean[indices[i-1]] - values[indices[i - k2 - 1]];
    mean[indices[i-1]] /= (n + k2 + 1 - i);
  }
  // normalize last entry
  mean[indices[n-1]] /= (k2 + 1);
}

// result should be allocated and zero
void dist(const float *p1, const float *p2, float *result) {
//  assert(p1 != NULL);
//  assert(p2 != NULL);
//  assert(result != NULL);
  for (int i = 0; i < dim; ++i) {
  //  printf("%f, %f\n", p1[i], p2[i]);
    *result += ((p1[i] - p2[i]) * (p1[i] - p2[i]));
  }
//  printf("%f\n", *result);
  *result = sqrt(fabs(*result));
}

//#define FIXED_NEIGHBOURS

/*
 * Global offset in dim dimensions.
 */
void calc_global_offset(int n, const float *xy, float *offset) {
 
  float mean[dim];
  for (int j = 0; j < dim; ++j) {
    mean[j] = avg(n, xy, j);
  }
  for (int j = 0; j < dim; ++j) {
    for (int i = 0; i < n; ++i) {
      offset[i*dim+j] = mean[j];
    }
  }
}

#define FIXED_WINDOW

#ifdef FIXED_WINDOW

/**
 * The number of neighbours is flexible.
 * The distance is fixed. 
 * Returns offset matched with corresponding xy values (not sorted).
 * Uses argsort / quicksort under the hood to calculate distances.
 */
void calc_offset(int n, const float *xy, float *offset) {

  // scalar distance per pair of points
  float *distances = (float*)calloc(sizeof(float),n*n);

  // calculate distances between all pairs of points (we later only need the points closer than window, so this can
  // be optimized)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      dist(&xy[i*dim], &xy[j*dim], &distances[j*n+i]);
    }
  }

  // for each point temporarily store indices to neighbours ordered from nearby to far away
  int *indices = (int*)malloc(sizeof(int)*n);

  // clear offsets (size n * dim)
  for (int i = 0; i < n*dim; ++i) {
    offset[i] = 0;
  }

  // fixed distance
  static float window = 2.0;

  for (int i = 0; i < n; ++i) {
    float *dist_for_i = &distances[i*n];
    argsort(n, dist_for_i, indices);

    int m_cnt = count_items_below_threshold(n, dist_for_i, indices, window);

    for (int j = 0; j < m_cnt; ++j) {
      offset[i*dim + 0] += xy[indices[j]*dim + 0]; 
      offset[i*dim + 1] += xy[indices[j]*dim + 1]; 
    }
    if (m_cnt > 0) {
      offset[i*dim + 0] /= (float)m_cnt;
      offset[i*dim + 1] /= (float)m_cnt;
    }
  }

  free(distances);
  distances = NULL;
}

#endif




#ifdef FIXED_NEIGHBOURS

/**
 * Currently, we sort for k/2 neighbours. However, the number of neighbours to be considered should be flexible.
 * The distance should be fixed and at around the maximum size of objects to be discovered. 
 */
// assumes offset is malloced
void calc_offset(int n, const float *xy, float *offset) {

  // scalar distance per pair of points
  float *distances = (float*)calloc(sizeof(float),n*n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      dist(&xy[i*dim], &xy[j*dim], &distances[j*n+i]);
    }
  }

  int *indices = (int*)malloc(sizeof(int)*n);

  // clear offsets (size n * dim)
  for (int i = 0; i < n*dim; ++i) {
    offset[i] = 0;
  }


  for (int i = 0; i < n; ++i) {
    float *dist_for_i = &distances[i*n];
    argsort(n, dist_for_i, indices);

    // number of neighbours
    static int k = n/2;

    for (int d = 0; d < dim; ++d) {
      for (int j = 0; j < k; ++j) {
        offset[i*dim + d] += xy[indices[j]*dim + d]; 
      }
      offset[i*dim + d] /= (float)k;
    }
  }

  free(distances);
  distances = NULL;
}
#endif

