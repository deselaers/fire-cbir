#include "diag.hpp"
#include "getpot.hpp"
#include "histogramfeature.hpp"
#include "sparsehistogramfeature.hpp"

int main(int argc, char **argv) {
  SparseHistogramFeature sh;
  sh.load(argv[1]);

  HistogramFeature h(sh);
  h.save(argv[2]);
}
