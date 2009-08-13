#include <iostream>
#include <limits>
#include "diag.hpp"
#include "jflib.hpp"
#include "gzstream.hpp"

using namespace std;

void USAGE() {
  cout << "USAGE: jf2arff <.jf file> <.arff file>";
  exit(1);
}

int main(int argc, char** argv) {

  if (argc != 3) {
    USAGE();
  }

  pair<DoubleVectorVector, IntVector> jfdata = readJF(argv[1]);
  int maxClass = 0;
  int minClass = numeric_limits<int>::max();
  if (jfdata.second.size() == 0) {
    ERR << ".jf file is empty or could not be read" << endl;
    exit(1);
  }
  for (uint i = 0; i < jfdata.second.size(); i++) {
    if (jfdata.second[i] < 0) {
      ERR << "did not expect a negative class label: " << jfdata.second[i] << "!" << endl;
      exit(1);
    }
    if (jfdata.second[i] > maxClass) {
      maxClass = jfdata.second[i];
    }
    if (jfdata.second[i] < minClass) {
      minClass = jfdata.second[i];
    }
  }
  if (minClass > 0) {
    ERR << "classes should be labeled starting from '0' !" << endl;
    ERR << "minimal class value found is '" << minClass << "'" << endl;
    exit(1);
  }

  ofstream os; os.open(argv[2], ios_base::out);
  if(!os.good()) {
    ERR << "Cannot write to " << argv[2] << endl;
    exit(20);
  }

  os << "@RELATION i6" << endl;
  os << "@ATTRIBUTE class {0";
  for (uint c = 1; c <= (uint) maxClass; c++) {
    os << "," << c;
  }
  os << "}" << endl;
  for (uint a = 0; a < jfdata.first[0]->size(); a++) {
    os << "@ATTRIBUTE attr" << a << " REAL" << endl;
  }
  os << "@DATA" << endl;

  for (uint i = 0; i < jfdata.first.size(); i++) {
    os << jfdata.second[i];
    for (uint j = 0; j < jfdata.first[i]->size(); j++) {
      os << "," << (*jfdata.first[i])[j];
    }
    os << endl;
  }

}
