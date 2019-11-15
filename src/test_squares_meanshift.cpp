#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <cassert>
#include <fstream>
#include <sys/stat.h> 
#include <emd_meanshift.h>

// Check the split from one to two clusters
//#define TWO_CLUSTERS

typedef struct _square_cfg {
  float scale[2];
  float translation[2];
} square_cfg_t;

void swap(const int i, const int j, float* pnts) {
  float tmp;
  for (int d = 0; d < dim; ++d) {
    tmp = pnts[j * dim + d];
    pnts[j * dim + d] = pnts[i * dim + d];
    pnts[i * dim + d] = tmp;
  }
}

void inverse(const int n, float *pnts) {
  for (int i = 0; i < n/2; ++i) {
    swap(i,n-1-i,pnts);
  }
}

void shuffle(const int n, float *pnts) {
  float tmpx, tmpy;
  for (int i = 0; i < n; ++i) {
    int j = rand() % n;
    swap(i,j,pnts);
  }
}

//#define DETERMINISTIC

void create_square(const int n, const square_cfg_t & square_cfg, float *pnts) {

  assert(dim == 2);

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
int main(int argc, char **argv) {

  // very simple way to check one argument "--method"
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << "--method METHOD" << std::endl;
    std::cout << "  wrong number of arguments" << std::endl;
    exit(1);
  }

  std::string method_key = std::string(argv[1]);
  std::string method_value = std::string(argv[2]);

  if (method_key != "--method") {
    std::cout << "Usage: " << argv[0] << "--method METHOD" << std::endl;
    std::cout << "  wrong second argument: " << method_key << std::endl;
    exit(1);
  }
  
  if (method_value == "default") {
    std::cout << "Run default emd algorithm" << std::endl;
  }
  else if (method_value == "shift") {
    std::cout << "Run shift-emd algorithm where k is determined in a flexible way." << std::endl;
  }
  else if (method_value == "global") {
    std::cout << "Run global shift-emd algorithm." << std::endl;
  } 
  else {
    std::cout << "Usage: " << argv[0] << "--method METHOD" << std::endl;
    std::cout << "  unknown method: " << method_value << std::endl;
    exit(1);
  }

  //int seed = 1234567;
  int seed = 11111;
 // seed = time(NULL);
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

#ifdef TWO_CLUSTERS
  std::string ofolder = "two_clusters";
#else
  std::string ofolder = "single_cluster";
#endif
  if (mkdir(ofolder.c_str(), 0777) == -1) {
    if( errno == EEXIST ) {
      std::cout << "Output directory exists, files in it might be overwritten." << std::endl;
    } else {
      std::cerr << "Could not create output directory: " << ofolder << std::endl;
      exit(2);
    }
  }
  std::string ofoldermethod = ofolder + '/' + method_value;
  if (mkdir(ofoldermethod.c_str(), 0777) == -1) {
    if( errno == EEXIST ) {
      std::cout << "Output directory exists, files in it might be overwritten." << std::endl;
    } else {
      std::cerr << "Could not create output directory: " << ofoldermethod << std::endl;
      exit(2);
    }
  }
  std::string fprefix = std::string(ofolder + '/' + method_value + '/');

  std::cout << "Write point cloud 1 to file" << std::endl;
  std::ofstream cloud1;
  cloud1.open(fprefix + "cloud1.txt");
  for (int i = 0; i < n; ++i) {
    cloud1 << std::fixed << std::setprecision(5) << xy1[i*dim] << "," << xy1[i*dim+1] << std::endl;
  }
  cloud1.close();

  std::cout << "Write point cloud 2 to file" << std::endl;
  std::ofstream cloud2;
  cloud2.open(fprefix + "cloud2.txt");
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

  if (method_value == "default") {
    emd_standard(b, n, m, xy1, xy2, match, offset1, offset2);
  }
  else if (method_value == "shift") {
    emd_meanshift(b, n, m, xy1, xy2, match, offset1, offset2);
  } 
  else if (method_value == "global") {
    emd_global_offset(b, n, m, xy1, xy2, match, offset1, offset2);
  }

  std::cout << "Write offset cloud 1 to file" << std::endl;
  std::ofstream foffset1;
  foffset1.open(fprefix + "offset1.txt");
  for (int i = 0; i < n; ++i) {
    foffset1 << std::fixed << std::setprecision(5) << offset1[i*dim] << "," << offset1[i*dim+1] << std::endl;
  }
  foffset1.close();

  std::cout << "Write offset cloud 2 to file" << std::endl;
  std::ofstream foffset2;
  foffset2.open(fprefix + "offset2.txt");
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
  fmatch.open(fprefix + "match.txt");
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

  if (method_value == "global") {
    emd_mean_costs_global_offset(b, n, m, xy1, xy2, match, offset1, offset2, cost);
  }
  else {
    emd_costs(b, n, m, xy1, xy2, match, offset1, offset2, cost);
  }
  std::cout << "Cost: " << cost[0]/n << std::endl;

  delete [] xy1;
  delete [] xy2;
  delete [] match;
  delete [] offset1;
  delete [] offset2;
  delete [] cost;

  return EXIT_SUCCESS;
}
