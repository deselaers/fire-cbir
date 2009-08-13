#include "sparsehistogramfeature.hpp"
#include "dist_jsd.hpp"
#include "string.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <sys/time.h>

using namespace std;

// for testing
int main(int argc, char** argv) {

  if (argc < 4) {
    ERR << "Usage: testsparsehistogramfeature <histofile 1> <histofile 2> <# iter> <smooth factor>" << endl;
    exit(1);
  }

  SparseHistogramFeature* shf1 = new SparseHistogramFeature();

  DBG(10) << "loading " << argv[1] << endl;
  shf1->load(argv[1]);
  DBG(10) << "loaded " << endl;

  SparseHistogramFeature* shf2 = new SparseHistogramFeature();
  DBG(10) << "loading " << argv[2] << endl;
  shf2->load(argv[2]);
  DBG(10) << "loaded " << endl;

  JSDDistance dist(atof(argv[4]));
  double sum = 0.0;
  double distance = 0.0;

  struct timeval startTime, endTime;
  gettimeofday(&startTime, NULL);

  int iter = atoi(argv[3]);
  for (int i = 0; i < iter; i++) {
    distance = dist.distance(shf1, shf2);
    sum += distance;
  }

  gettimeofday(&endTime, NULL);
  long millisecs = (endTime.tv_sec - startTime.tv_sec) * 1000;
  if (endTime.tv_usec >= startTime.tv_usec) {
    millisecs += (endTime.tv_usec - startTime.tv_usec) / 1000;
  } else {
    millisecs += (startTime.tv_usec - endTime.tv_usec) / 1000;
    millisecs -= 1000;
  }

	delete shf1;
	delete shf2;

  DBG(10) << "distance is " << distance << endl;
  DBG(10) << "calculation took " << millisecs << " ms." << endl;
  return (int) sum;

}
