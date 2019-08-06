#pragma once

// currently, experiment with objects in 2D
static const int dim = 2;

/**
 * Adaptation of EMD algorithm that calculates offsets and takes these into account.
 */
void emd_meanshift(int b, int n, int m, const float * xy1, const float * xy2, float * match, 
    float * offset1, float * offset2);

/**
 * Paper: https://arxiv.org/pdf/1612.00603.pdf
 *
 * d_EMD(S1,S2) = min_\phi \sum_{x \in S1} || x - phi(x) ||_2
 *
 * The earth mover distance calculates the Euclidean distance || x - phi(x) ||_2. The optimal bijection \phi finds
 * the closest point in S2 with respect to the given point in S1. This is called the assignment problem.
 *
 * Here a (1 + \epsilon) approximation scheme is used for the assignment problem, also called the bipartite perfect
 * matching problem.
 *
 * Paper: https://ieeexplore.ieee.org/abstract/document/4048607/ (Bertsekas, 1985)
 *
 * This does not seem to be the case... The implementation does not look like an auction. It has weights, it might
 * be something like this:  http://www.columbia.edu/~cs2035/courses/ieor8100.F18/GabTar.pdf. It is not clear which
 * implementation has been used for the matching.
 *
 * Approximate the match using some kind of Earth Mover's Distance / 1-Wasserstein Distance. 
 *
 * We find the matching point for each element in xy1 in the matrix xy2.
 *
 * Output: match matrix of size b x n x m.
 *         offset is twice the size of the match, for every pair of points it indicates how much to shift the 
 *         second point unto the first before calculating their distance
 *
 * @param b        number of batches
 * @param n        number of points in point cloud 1 (batch)
 * @param m        number of points in point cloud 2 (batch)
 * @param xy1      the xy coordinates in point cloud 1 in format [x0 y0 z0 x1 y1 z1 ... xn yn zn]
 * @param xy2      the xy coordinates in point cloud 2 in format [x0 y0 z0 x1 y1 z1 ... xn yn zn]
 * @param match    result, zero matrix with positive values for transportation between points in xy1 and xy2
 * @param offset1  matrix of size n * dim. 
 * @param offset2  matrix of size m * dim. 
 */
void approxmatch_cpu(int b, int n, int m, const float * xy1, const float * xy2, float * match, 
    float * offset1, float * offset2);

/**
 * The cost function. We calculate the cost for each item in the batch b.
 * Input xyz1 is of dimension b x n.
 * Input xyz2 is of dimension b x m.
 * Input match is of dimension b x n x m. It is 1 if the points in xyz1 and xyz2 match.
 *
 * For each matching point we calculate the Euclidian distance. Note that this is the 1-Wasserstein distance. 
 * The distance metric is Euclidean and it is not squared p=2 or cubed p=3, or otherwise.
 * The cost is just the sum of Euclidean distances.
 *
 * If b = 1, n is total number of points in point cloud 1. If b = 2, n should be half of that.
 *
 * @param b        number of batches
 * @param n        number of points in point cloud 1 (batch)
 * @param m        number of points in point cloud 2 (batch)
 * @param xy1      the xy coordinates in point cloud 1 in format [x0 y0 x1 y1 ... xn yn zn]
 * @param xy2      the xy coordinates in point cloud 2 in format [x0 y0 x1 y1 ... xn yn zn]
 * @param match    zero matrix with only 1s when points in xyz1 and xyz2 match
 * @param offset1  matrix of size n * dim. 
 * @param offset2  matrix of size m * dim. 
 * @param cost     result, for each matching point, calculate euclidean distance and calculate the overall sum
 */
void matchcost_cpu(int b, int n, int m, float * xy1, float * xy2, float * match,
    float * offset1, float * offset2, float * cost);
