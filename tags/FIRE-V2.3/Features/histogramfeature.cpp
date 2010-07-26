/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

FIRE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <sstream>
#include <string>
#include <vector>
#include "gzstream.hpp"
#include "diag.hpp"
#include "histogramfeature.hpp"

using namespace std;


HistogramFeature::HistogramFeature(uint counter, const vector<uint> & bins) : VectorFeature(bins.size()), bins_(bins),steps_(1,bins.size()),min_(1),max_(1),stepSize_(1),size_(bins.size()),dim_(1),counter_(counter){
  for(uint i=0;i<bins_.size();++i) {
    data_[i]=double(bins_[i])/double(counter_);
  }
}
  


HistogramFeature::HistogramFeature(uint size) : VectorFeature(size), bins_(size),steps_(1,size),min_(1),max_(1),stepSize_(1),size_(size),dim_(1),counter_(0)  { type_=FT_HISTO; }

HistogramFeature::HistogramFeature(::std::vector<uint> steps) : steps_(steps), 
                                                                min_(steps.size(),0.0), 
                                                                max_(steps.size(),1.0),
                                                                stepSize_(steps.size()),
                                                                dim_(steps.size()),
                                                                counter_(0) {
  type_=FT_HISTO;
  size_=1;
  for(uint i=0;i<steps_.size();++i) {
    size_*=steps_[i];
  }
  bins_=vector<uint>(size_);
  data_=vector<double>(size_);
}


HistogramFeature::HistogramFeature(SparseHistogramFeature sh) : steps_(sh.steps()), 
                                                                min_(sh.min()), 
                                                                max_(sh.max()), 
                                                                dim_(sh.steps().size()),
                                                                counter_(sh.counter())
                                                                
{
  type_=FT_HISTO;
  initStepsize();
  size_=1;
  for(uint i=0;i<steps_.size();++i) {
    size_*=steps_[i];
  }

  vector<double> tmp; sh.getAllBinValues(tmp);
  
  data_=vector<double>(size_);
  bins_=vector<uint>(size_);

  for(uint i=0;i<tmp.size();++i) {
    data_[i]=tmp[i];
    bins_[i]=int(tmp[i]*counter_);
  }
}

HistogramFeature::HistogramFeature() : VectorFeature(),bins_(0),steps_(1),min_(1),max_(1),stepSize_(1),size_(0),dim_(1),counter_(0) {type_=FT_HISTO; }

HistogramFeature::~HistogramFeature() {
  bins_.clear();
  data_.clear();
  steps_.clear();
}


void HistogramFeature::initStepsize() {
  stepSize_=vector<double>(dim_);
  for(uint i=0;i<dim_;++i) {
    stepSize_[i]=(max_[i]-min_[i])/double(steps_[i]);
  }
}

uint HistogramFeature::posToBin(const vector<uint>& pos) const {
  uint p=0;
  if(pos.size()!=steps_.size()) {
    ERR << "Invalid position: not right dimensionality." << endl;
    p=bins_.size()+1;
  } else {
    uint m=1;
    for(uint i=0;i<pos.size();++i) {
      if(pos[i]<steps_[i]) {
        p=p+pos[i]*m;
        m*=steps_[i];
      } else {
        ERR << "Invalid position: out of range." << endl;
        p=bins_.size()+1;
      }
    }
  }
  return p;
}

double HistogramFeature::operator()(const ::std::vector<uint>& pos) const {
  uint p=this->posToBin(pos);
  if(p < bins_.size()) {
    return data_[p];
  } else {
    return -1.0;
  }
}

  
const uint HistogramFeature::size() const{
  return bins_.size();
}

const unsigned long int HistogramFeature::calcBinarySize() const{
  unsigned long int bsize = 0;
  // the numbers for dimension and counter
  bsize += 2*sizeof(uint);
  // add up the size of the vectors
  bsize += (unsigned long int)steps_.size() * (unsigned long int)sizeof(uint);
  bsize += (unsigned long int)min_.size() * (unsigned long int)sizeof(double);
  bsize += (unsigned long int)max_.size() * (unsigned long int)sizeof(double);
  bsize += (unsigned long int)bins_.size() * (long unsigned int)sizeof(uint);
  return bsize;
}

const uint& HistogramFeature::bin(const uint idx) const{
  return bins_[idx];
}

void HistogramFeature::feedbin(const uint idx) {
  ++counter_;
  ++bins_[idx];
  data_[idx]=double(bins_[idx])/double(counter_);
}

const uint& HistogramFeature::bin(const ::std::vector<uint>& pos) const {
  return bins_[posToBin(pos)];
}

void HistogramFeature::feedbin(const ::std::vector<uint>& pos) {
  ++counter_;
  uint p=posToBin(pos);
  ++bins_[p];
  data_[p]=double(bins_[p])/double(counter_);
}



vector<uint> HistogramFeature::pointToPos(const vector<double>& point) const {
  vector<uint> pos(dim_,0);
  
  BLINK(55) << "position: " ;  
  for(uint i=0;i<dim_;++i) {
    pos[i]=int(point[i]/stepSize_[i]);
    if(pos[i]>=steps_[i]) pos[i]=steps_[i]-1;
    //BLINK(55) << point[i] << " -> " << pos[i] << "\t ";
    //cout << point[i] << " -> " << pos[i] << endl;
  }
  BLINK(55) << endl;
  
  return pos;
}

void HistogramFeature::feed(const ::std::vector<double>& inData) {
  ::std::vector<uint> pos = pointToPos(inData);
  //cout << "binned to " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  this->feedbin(pos);
}

//load old fire v0.1 version histograms
void HistogramFeature::loadOld(const ::std::string &filename) {
  dim_=1;
  igzstream is; is.open(filename.c_str());
  if(!is.good()) {
    ERR << "Unable to open histogram (old style)'"<< filename << "'."<< endl;
    return;
  } else {
    string line;

    while(!is.eof()) {
      if('#'==is.peek()) { // comment line
        getline(is,line); 
      } else {
        getline(is,line);
        if(!is.eof()) {
          istringstream iss(line);
          string token;
          iss >> token;
          if(token=="dim") {
            iss >> dim_;
          } else if(token=="steps") {
            uint tmp;
            iss >> tmp;
            steps_=vector<uint>(dim_,tmp);
            size_=1;
            for(uint i=0;i<steps_.size();++i) {
              size_*=steps_[i];
            }
          } else if(token=="min") {
            min_=vector<double>(0);
            for(uint i=0;i<dim_;++i) {
              double tmp;
              iss >> tmp;
              min_.push_back(tmp);
            }
            if(dim_!=min_.size()) {
              ERR << "min "<<min_.size() << "  and dim " << dim_ << " inconsistent while reading: " <<filename << endl;
            }
          } else if(token=="max") {
            max_=vector<double>(0);
            for(uint i=0;i<dim_;++i) {
              double tmp;
              iss >> tmp;
              max_.push_back(tmp);
            }
            if(dim_!=max_.size()) {
              ERR << "max " << max_.size() << " and dim " << dim_ << " inconsistent while reading: "<< filename << endl;
            }
          } else if (token=="stepsize") {
            //forget it
          } else if(token=="nOfBins") {
            uint test;
            iss >> test;
            if(test != size_) {
              ERR << "nofBins " << test << " and size " << size_ << " inconsistent while reading: " <<filename<< endl;
            }
          } else if (token=="counter") {
            iss >> counter_;
          } else if(token=="bins") {
            DBG(30) << "OK, now reading bins" << endl;
            bins_=vector<uint>(size_);
            data_=vector<double>(size_);
            for(unsigned int i=0;i<size_;++i) {
              uint number, content;
              is >> number >> content;
              if(number != i) {
                ERR << "While reading bins: number " << number << " != i " << i << " while reading: " << filename << endl;
              }
              bins_[i]=content;
              data_[i]=double(content)/double(counter_);
            }
          } else if (token!="") {
            ERR << "Unknown token in line: " << line << " while reading " << filename <<endl;
          } 
          iss.clear();
        }
      }
    }
  }
  is.close();
}

bool HistogramFeature::read(istream &is) {
  istringstream iss;
  string line,token;

  // MagickNumber
  getline(is,line);
  if(line!="FIRE_histogram") {
    ERR << "Not reading a valid Histogram:" << endl;
    return false;
  }

  // Comments
  while('#'==is.peek()) { // comment lines
    getline(is,line); 
  }
  
  // dim
  token=""; getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="dim") {
    iss >> dim_;
  } else {
    ERR << "Expected 'dim', got '" << line <<"'." << endl;
    return false;
  }

  // counter
  token=""; getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="counter"){
    iss >> counter_;
  } else {
    ERR << "Expected 'counter', got '" << token <<"'." << endl;
    return false;
  }

  if (counter_==0) { 
    counter_=1;
    ERR << "Histogram without data points. Setting counter_:=1 -> histogram has zero in all bins." << endl;
    return false;
  }
  
  // steps
  getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="steps") {
    steps_=vector<uint>(0);
    size_=1;
    while(!iss.eof()) {
      uint tmp;
      iss >> tmp;
      steps_.push_back(tmp);
      size_*=tmp;
    }
    if(dim_!=steps_.size()) {
      // commented out due to bug in histogramfeature
      //ERR << "steps " << steps_.size() << " and dim " << dim_ << " inconsistent." << endl;
    }
  } else {
    ERR << "Expected 'steps', got '" << line <<"'." << endl;
    return false;
  }

  // min
  getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="min") {
    min_=vector<double>(0);
    while(!iss.eof()) {
      double tmp;
      iss >> tmp;
      min_.push_back(tmp);
    }
    if(dim_!=min_.size()) {
      ERR << "min "<<min_.size() << " and dim " << dim_ << " inconsistent." << endl;
    }
  } else {
    ERR << "Expected 'min', got '" << line <<"'." << endl;
    return false;
  }

  //max
  getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="max") {
    max_=vector<double>(0);
    while(!iss.eof()) {
      double tmp;
      iss >> tmp;
      max_.push_back(tmp);
    }
    if(dim_!=max_.size()) {
      ERR << "max " << max_.size() << " and dim " << dim_ << " inconsistent." << endl;
    }
  } else {
    ERR << "Expected 'max', got '" << line <<"'." << endl;
    return false;
  }
  
  // data
  getline(is,line); iss.clear(); iss.str(line); iss >> token;
  if(token=="data"){
    bins_=vector<uint>(0);
    data_=vector<double>(0);
    while(!iss.eof()) {
      uint tmp;
      iss >> tmp;
      bins_.push_back(tmp);
      data_.push_back(double(tmp)/double(counter_));
    }
    if(size_!=bins_.size()) {
      // commented out due to bug
      //ERR << "data " << bins_.size() << " and size " << size_ << " inconsistent." << endl;
    }
  } else {
    ERR << "Expected 'data', got '" << line <<"'." << endl;
    return false;
  }
  this->initStepsize();
  iss.clear();
  return true;
}

bool HistogramFeature::readBinary(istream &is){
	if(!is.good()){
		return false;
	}	
	// read dimension marker
	is.read((char*)&dim_,sizeof(uint));
	if(is.fail() || is.eof()){
		return false;
	}
	// read counter
	is.read((char*)&counter_,sizeof(uint));
	if(is.fail() || is.eof()){
		return false;
	}
	if (counter_==0) { 
      counter_=1;
      ERR << "Histogram without data points. Setting counter_:=1 -> histogram has zero in all bins." << endl;
    }
    
    // read steps
    // note that due to binary restrictions the amount of steps must be
    // consistent with the dimension of the histogram
    steps_ = vector<uint>(0);
    size_=1;
    for(uint i = 0; i< dim_;++i){
    	uint tmp;
    	is.read((char*)&tmp,sizeof(uint));
    	if(is.fail() || is.eof()){
    		return false;
    	}
    	steps_.push_back(tmp);
    	size_*=tmp;            
   	 }
	
    
    // read min
	min_=vector<double>(0);
	for (uint i = 0; i < dim_; ++i) {
		double tmp;
		is.read((char*)&tmp,sizeof(double));
		if(is.fail() || is.eof()){
			return false;
		}
		min_.push_back(tmp);
	}
        
    // read max
    max_=vector<double>(0);
    for (uint i = 0; i < dim_; ++i) {
    	double tmp;
     	is.read((char*)&tmp,sizeof(double));
     	if(is.fail() || is.eof()){
     		return false;
     	}
      	max_.push_back(tmp);
    }
        
    // read data
	bins_=vector<uint>(0);
	data_=vector<double>(0);
	for(uint i = 0; i< size_;++i) {
	   uint tmp;
	   is.read((char*)&tmp,sizeof(uint));
	   if(is.fail() || (i<size_-1 && is.eof())){
	   		return false;
	   }
	   bins_.push_back(tmp);
	   data_.push_back(double(tmp)/double(counter_));
	}
	this->initStepsize();
	return true;
}

void HistogramFeature::write(ostream &os) {
  os << "FIRE_histogram" << endl
     << "# Histogram file saved for Fire V2" << endl
     << "dim " << dim_ << endl
     << "counter " << counter_ << endl
     << "steps";
  for(uint i=0;i<steps_.size();++i) {os << " "<< steps_[i];} os << endl;
  os << "min";
  for(uint i=0;i<min_.size();++i) {os << " "<< min_[i];} os << endl;
  os << "max";
  for(uint i=0;i<max_.size();++i) { os << " "<< max_[i]; } os << endl;
  os << "data";
  for(uint i=0;i<bins_.size();++i) { os << " " << bins_[i]; } os << endl;
}

void HistogramFeature::writeBinary(ostream &os) {
	os.write((char*)&dim_,sizeof(uint));
	os.write((char*)&counter_,sizeof(uint));
	for(uint i=0; i< steps_.size();++i) {os.write((char*)&steps_[i],sizeof(uint)); }
	for(uint i=0; i< min_.size();++i) {os.write((char*)&min_[i],sizeof(double)); }
	for(uint i=0; i< max_.size();++i) {os.write((char*)&max_[i],sizeof(double)); }
	for(uint i=0; i< bins_.size();++i) {os.write((char*)&bins_[i],sizeof(uint)); }
}
