#include <string>
#include "gzstream.hpp"
#include "basefeature.hpp"

using namespace std;

bool BaseFeature::load(const ::std::string &filename) {
  DBG(30) << "Loading from filename '" << filename << "'." << endl;
  igzstream is; // if igzstream is constructed with the file to be
                // opened, the state is wrong! Thus we have
                // construction and opening in two different lines
  is.open(filename.c_str());
  if(!is.good()) {
    ERR << "Canot open '" << filename << " for reading." << endl;
    return false;
  } else {
    if(not this->read(is)) {
      ERR << "Problem when reading '" << filename << "'." << endl;
      return false;
    }
  }
  DBG(40) << "Loaded from filename '" << filename << "'." << endl;
  is.close();
  return true;
}

void BaseFeature::save(const ::std::string &filename) {
  DBG(30) << "Writing to '" << filename << "'." << endl;
  ogzstream os; // if ogzstream is constructed with the file to be
                // opened, the state is wrong! Thus we have
                // construction and opening in two different lines
  os.open(filename.c_str());
  if(!os.good()) {
    ERR << "Cannot open '" << filename << "' for writing." << endl;
  } else {
    this->write(os);
  }
  DBG(40) << "Saved to '" << filename << "'." << endl;
}
