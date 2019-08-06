#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <cassert>
#include <fstream>

#include <emd_meanshift.h>

extern "C" {
#include <sort_indices.h>
}

void emd_meanshift(int b, int n, int m, const float * xy1, const float * xy2, float * match, 
    float * offset1, float * offset2) {

  // here offset is just per point, not per pair
  calc_offset(n, xy1, offset1);
  calc_offset(m, xy2, offset2);

  approxmatch_cpu(b, n, m, xy1, xy2, match, offset1, offset2);
}

void approxmatch_cpu(int b, int n, int m, const float * xy1, const float * xy2, float * match, 
    float * offset1, float * offset2) {

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

void matchcost_cpu(int b, int n, int m, float * xy1, float * xy2, float * match, float * offset1, float * offset2, 
    float * cost) {

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
