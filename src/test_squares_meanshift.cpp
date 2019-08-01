#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <cassert>
#include <fstream>

extern "C" {
#include <sort_indices.h>
}

// play around with 2D rather than 3D
static const int dim = 2;

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
 * @param offset   matrix of size m * n * dim. 
 *
 * TODO: change to offset1 and offset2
 * 
 */
void approxmatch_cpu(int b,int n,int m,const float * xy1,const float * xy2,float * match, float * offset1, float *offset2) {
  
  // here offset is just per point, not per pair
  calc_offset(n, xy1, offset1);
  calc_offset(m, xy2, offset2);

  for (int i=0;i<b;i++){

    // decompose in such way that one factor is 1 and the other factor defines how often the cloud point fits
    int factorl=std::max(n,m)/n;
    int factorr=std::max(n,m)/m;

    // saturation says something about convergence, initialize at factor L and factor R.
    std::vector<double> saturatedl(n,double(factorl)), saturatedr(m,double(factorr));
    // weights for each pair of points
    std::vector<double> weight(n*m);
    // init match matrix to 0
    for (int j=0;j<n*m;j++)
      match[j]=0;
    // iterate over 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2
    for (int j=8;j>=-2;j--){
      // level is then -65536, -16384, ..., -4, -1, -1/4, -1/16, the latter of which is set to 0
      double level=-powf(4.0,j);
      if (j==-2)
	level=0;
      for (int k=0;k<n;k++){
	double x1=xy1[k*dim+0];
	double y1=xy1[k*dim+1];
	// this iterates over all points, they all have the same 
	for (int l=0;l<m;l++){
	  double dx=(x1-xy2[l*dim+0]) - (offset1[k*dim+0] - offset2[l*dim+0]);
	  double dy=(y1-xy2[l*dim+1]) - (offset1[k*dim+1] - offset2[l*dim+1]);

	  // this is not sparse, that's why almost any point wants to contribute to other points
	  // even if there distance is not small
	  double dist = dx*dx+dy*dy;
	  //					double dist = abs(dx+dy);
	  //					dist = sqrtf(dist);
	  weight[k*m+l]=expf(level*dist)*saturatedr[l];
	}
      }
      // vector ss is sum for each l
      std::vector<double> ss(m,1e-9);
      for (int k=0;k<n;k++){
	double s=1e-9;
	// sum all weights
	for (int l=0;l<m;l++){
	  s+=weight[k*m+l];
	}
	// normalize with sum and multiply each point in k with saturation L
	for (int l=0;l<m;l++){
	  weight[k*m+l]=weight[k*m+l]/s*saturatedl[k];
	}
	// sum again for each point in l
	for (int l=0;l<m;l++) {
	  ss[l]+=weight[k*m+l];
	}
      }
      // normalize now over l
      for (int l=0;l<m;l++){
	double s=ss[l];
	double r=std::min(saturatedr[l]/s,1.0);
	ss[l]=r;
      }

      // bidding phase...
      // we should somehow find the biggest value v_i and then the second best entry w_i
      // vector ss2 is yet another sum
      std::vector<double> ss2(m,0);
      for (int k=0;k<n;k++){
	double s=0;
	for (int l=0;l<m;l++){
	  // we multiply the weights with ss
	  weight[k*m+l]*=ss[l];
	  // we add them to the sum s
	  s+=weight[k*m+l];
	  // we add them also to the sum ss2
	  ss2[l]+=weight[k*m+l];
	}
	// here we calculate saturated L as saturated L minus s
	saturatedl[k]=std::max(saturatedl[k]-s,0.0);
      }
      // write match matrix by adding weight, how is it only 0 or 1, it does not seem so.
      for (int kk=0;kk<n*m;kk++) {
	match[kk]+=weight[kk];
      }

      // saturation of R minus ss2
      for (int l=0;l<m;l++){
	saturatedr[l]=std::max(saturatedr[l]-ss2[l],0.0);
      }
    }
    xy1+=n*dim;
    xy2+=m*dim;
    match+=n*m;
  }
}

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
 * @param cost     result, for each matching point, calculate euclidean distance and calculate the overall sum
 */
void matchcost_cpu(int b,int n,int m,float * xy1,float * xy2,float * match,float *offset1, float *offset2, float * cost){
	for (int i=0;i<b;i++){
		double s=0;
		for (int j=0;j<n;j++)
			for (int k=0;k<m;k++){
				float x1=xy1[j*dim+0];
				float y1=xy1[j*dim+1];
				//float dx=x1 - xy2[k*dim+0] - offset[(j*m+k)*dim+0];
				//float dy=y1 - xy2[k*dim+1] - offset[(j*m+k)*dim+1];
				float dx=(x1 - xy2[k*dim+0]) - (offset1[j*dim+0] - offset2[k*dim+0]);
				float dy=(y1 - xy2[k*dim+1]) - (offset1[j*dim+1] - offset2[k*dim+1]);
				float d=sqrtf(dx*dx+dy*dy)*match[j*m+k];
				s+=d;
			}
		cost[0]=s;
		xy1+=n*dim;
		xy2+=m*dim;
		match+=n*m;
		cost+=1;
	}
}

typedef struct _square_cfg {
	float scale[2];
	float translation[2];
} square_cfg_t;

void swap(const int i, const int j, float* pnts) {
	float tmpx, tmpy;
	tmpx = pnts[j * dim + 0];
	tmpy = pnts[j * dim + 1];
	pnts[j * dim + 0] = pnts[i * dim + 0];
	pnts[j * dim + 1] = pnts[i * dim + 1];
	pnts[i * dim + 0] = tmpx;
	pnts[i * dim + 1] = tmpy;

}

void inverse(const int n, float *pnts) {
	for (int i = 0; i < n/2; ++i) {
		swap(i,n-1-i,pnts);
	}
}

void shuffle(const int n, float *pnts) {
	int dim = 2;
	float tmpx, tmpy;
	for (int i = 0; i < n; ++i) {
	 	int j = rand() % n;
		swap(i,j,pnts);
//		std::cout << "swap " << i << " and " << j << std::endl;
	}
}

//#define DETERMINISTIC

void create_square(const int n, const square_cfg_t & square_cfg, float *pnts) {
	
	int dim = 2;

	float *squareRnd = new float[n];

	for (int i = 0; i < n; ++i) {
	 	squareRnd[i] = drand48() - 0.5;
	 	
//	 	or make deterministic, sq = [-0.5 ... 0.5 -0.5 ... 0.5]
//	 	counts corners double
#ifdef DETERMINISTIC
		squareRnd[i] = (i%(n/4)*1.0 / (n/4-1)) - 0.5;
#endif
//		squareRnd[i] = (i%(n)*1.0 / (n-1)) - 0.5;
	}

	/* // single line
	for (int i = 0; i < n; ++i) {
		pnts[i * dim + 0] = squareRnd[i] * square_cfg.scale[0] + square_cfg.translation[0];
		pnts[i * dim + 1] =         -0.5 * square_cfg.scale[1] + square_cfg.translation[1];
	} */

	for (int i = 0; i < n/4; ++i) {
		pnts[i * dim + 0] = squareRnd[i] * square_cfg.scale[0] + square_cfg.translation[0];
		pnts[i * dim + 1] =         -0.5 * square_cfg.scale[1] + square_cfg.translation[1];
	}
	for (int i = n/4; i < n/2; ++i) {
		pnts[i * dim + 0] =         -0.5 * square_cfg.scale[0] + square_cfg.translation[0];
		pnts[i * dim + 1] = squareRnd[i] * square_cfg.scale[1] + square_cfg.translation[1];
	}
	for (int i = n/2; i < (n*3)/4; ++i) {
		pnts[i * dim + 0] = squareRnd[i] * square_cfg.scale[0] + square_cfg.translation[0];
		pnts[i * dim + 1] =          0.5 * square_cfg.scale[1] + square_cfg.translation[1];
	}
	for (int i = (n*3)/4; i < n; ++i) {
		pnts[i * dim + 0] =          0.5 * square_cfg.scale[0] + square_cfg.translation[0];
		pnts[i * dim + 1] = squareRnd[i] * square_cfg.scale[1] + square_cfg.translation[1];
	}

#ifndef DETERMINISTIC
	shuffle(n, pnts);
#endif

	delete [] squareRnd;
}

/**
 * Test approximation of earth mover's distance
 */
int main() {
	
  printf("Algorithm where k is determined in a flexible way.\n");

	int seed = 1213;
	seed = time(NULL);
	srand48(seed);

	float *xy1, *xy2, *match, *cost, *offset1, *offset2;
	int b = 1;
	// fix number of points to be the same number in cloud 1 and 2
	int n = 4*20;
	// for now, only allow multiple of 4
	assert(n % 4 == 0);
	int m = n;
	xy1 = new float[n*dim];
	xy2 = new float[m*dim];
	match = new float[m*n];
//	offset = new float[m*n*dim];
	offset1 = new float[n*dim];
	offset2 = new float[m*dim];

	cost = new float[1]; // batch of size 1

	square_cfg_t squareA_cfg;
	squareA_cfg.scale[0] = 1;
	squareA_cfg.scale[1] = 1;
	squareA_cfg.translation[0] = 0;
	squareA_cfg.translation[1] = 0;
	
	create_square(n, squareA_cfg, xy1);
	
	square_cfg_t squareB_cfg;
	squareB_cfg.scale[0] = 1;
	squareB_cfg.scale[1] = 1;
	squareB_cfg.translation[0] = 3;
	squareB_cfg.translation[1] = 2;

#define TWO_CLUSTERS

#ifdef TWO_CLUSTERS
	square_cfg_t squareC_cfg;
	squareC_cfg.scale[0] = 1;
	squareC_cfg.scale[1] = 1;
	squareC_cfg.translation[0] = -4;
	squareC_cfg.translation[1] = -5;

	create_square(n/2, squareB_cfg, xy2);
	create_square(n/2, squareC_cfg, xy2+n);
#else
	create_square(n, squareB_cfg, xy2);
#endif

	std::cout << "Write point cloud 1 to file" << std::endl;
	std::ofstream cloud1;
	cloud1.open("cloud1.txt");
	for (int i = 0; i < n; ++i) {
		cloud1 << std::fixed << std::setprecision(5) << xy1[i*dim] << "," << xy1[i*dim+1] << std::endl;
	}
	cloud1.close();

	std::cout << "Write point cloud 2 to file" << std::endl;
	std::ofstream cloud2;
	cloud2.open("cloud2.txt");
	for (int i = 0; i < n; ++i) {
		cloud2 << std::fixed << std::setprecision(5) << xy2[i*dim] << "," << xy2[i*dim+1] << std::endl;
	}
	cloud2.close();

	std::cout << std::endl;
	//std::cout << "Distances:" << std::endl;
	for (int k=0;k<n;k++){
		double x1=xy1[k*dim+0];
		double y1=xy1[k*dim+1];
		for (int l=0;l<m;l++){
			double x2=xy2[l*dim+0];
			double y2=xy2[l*dim+1];
			double dist = ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	//		std::cout << std::fixed << std::setw(4) << std::setfill('0') << std::setprecision(1) << dist << ' ' ;
		}
	//	std::cout << std::endl;
	}

	approxmatch_cpu(b, n, m, xy1, xy2, match, offset1, offset2);
	
	std::cout << "Write offset cloud 1 to file" << std::endl;
	std::ofstream foffset1;
	foffset1.open("offset1.txt");
	for (int i = 0; i < n; ++i) {
	  foffset1 << std::fixed << std::setprecision(5) << offset1[i*dim] << "," << offset1[i*dim+1] << std::endl;
	}
	foffset1.close();
	
	std::cout << "Write offset cloud 2 to file" << std::endl;
	std::ofstream foffset2;
	foffset2.open("offset2.txt");
	for (int i = 0; i < n; ++i) {
	  foffset2 << std::fixed << std::setprecision(5) << offset2[i*dim] << "," << offset2[i*dim+1] << std::endl;
	}
	foffset2.close();

	/*
	// but match should be something like
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			match[i*m+j] = (i == j);
		}
	}*/


	//std::cout << std::endl;
	//std::cout << "Match: " << std::endl;
	std::ofstream fmatch;
	fmatch.open("match.txt");
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
	//		std::cout << std::fixed << std::setprecision(4) << match[i*m+j] << ' ';
			fmatch << std::fixed << std::setprecision(4) << match[i*m+j] << ' ';

		}
	//	std::cout << std::endl;
		fmatch << std::endl;
	}
	//std::cout << std::endl;
	fmatch.close();

	int i = 0; int j = 0;
	std::cout << "Compare points with indices " << i << " and " << j << " in cloud 1 and 2" << std::endl;
	std::cout << "Cloud 1 p[x,y]=[" << xy1[i*dim] << "," << xy1[i*dim+1] << "]" << std::endl;
	std::cout << "Cloud 2 p[x,y]=[" << xy2[j*dim] << "," << xy2[j*dim+1] << "]" << std::endl;
	std::cout << "Match [i,j]=[" << i << "," << j << "]:" << match[i*m+j] << std::endl;
	std::cout << "Offset cloud 1 [i,j]=[" << i << "," << j << "] -> [" << offset1[j*dim+0] << "," << offset1[j*dim+1] << "]" << std::endl;
	std::cout << "Offset cloud 2 [i,j]=[" << i << "," << j << "] -> [" << offset2[j*dim+0] << "," << offset2[j*dim+1] << "]" << std::endl;
	i = 2; j = 2;
	std::cout << "Compare points with indices " << i << " and " << j << " in cloud 1 and 2" << std::endl;
	std::cout << "Cloud 1 p[x,y]=[" << xy1[i*dim] << "," << xy1[i*dim+1] << "]" << std::endl;
	std::cout << "Cloud 2 p[x,y]=[" << xy2[j*dim] << "," << xy2[j*dim+1] << "]" << std::endl;
	std::cout << "Match [i,j]=[" << i << "," << j << "]:" << match[i*m+j] << std::endl;
	std::cout << "Offset cloud 1 [i,j]=[" << i << "," << j << "] -> [" << offset1[j*dim+0] << "," << offset1[j*dim+1] << "]" << std::endl;
	std::cout << "Offset cloud 2 [i,j]=[" << i << "," << j << "] -> [" << offset2[j*dim+0] << "," << offset2[j*dim+1] << "]" << std::endl;
	
	
	for (int j = 0; j < m; ++j) {
//	  std::cout << "Cloud 2 p[x,y]=[" << xy2[j*dim] << "," << xy2[j*dim+1] << "]" << std::endl;
//	  std::cout << "Offset cloud 2 [i,j]=[" << i << "," << j << "] -> [" << offset2[j*dim+0] << "," << offset2[j*dim+1] << "]" << std::endl;
	}

	matchcost_cpu(b, n, m, xy1, xy2, match, offset1, offset2, cost);
	std::cout << "Cost: " << cost[0]/n << std::endl;

	delete [] xy1;
	delete [] xy2;
	delete [] match;
	delete [] offset1;
	delete [] offset2;
	delete [] cost;

	return EXIT_SUCCESS;
}
