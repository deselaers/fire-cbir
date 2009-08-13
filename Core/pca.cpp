#include <vector>
#include <fstream>
#include "diag.hpp"
#include "pca.hpp"
#include "gzstream.hpp"
#include "svd.hpp"
using namespace std;

PCA::PCA(int dim) {
  dim_=dim;
  mean_=vector<double>(dim_,0.0);
  covariance_=vector< vector<double> >(dim_,vector<double>(dim_,0.0));
  counter_=0;
}

void PCA::putData(const vector<double> &in) {
  ++counter_;
  double tmp;

  for(int i=0;i<dim_;++i) {
    mean_[i]+=in[i];
  }
  
  for(int i=0;i<dim_;++i) {
    tmp=in[i];
    for(int j=0;j<dim_;++j) {
      covariance_[i][j]+=tmp*in[j];
    }
  }
}
void PCA::save(const string &filename) const{
  ogzstream ofs; ofs.open(filename.c_str());
  if(!ofs.good()) {
    ERR << "[PCA::save] Cannot write PCA to '" << filename << "'."<< endl;
  } else {
    ofs << dim_ << " " <<dim_ << endl;
    ofs << "mean "  ;
    for(int i=0;i<dim_;++i) {
      ofs << mean_[i] << " "; }
    ofs << endl;
    
    for(int y=0;y<dim_;++y) {
      ofs << y;
      for(int x=0;x<dim_;++x) {
        ofs << " " << covariance_[y][x];
      }
      ofs << endl;
    }
    ofs << -1 << endl;
  }
}

bool PCA::load(const string &filename) {
  igzstream ifs; ifs.open(filename.c_str());
  int dimy, dimx, posy;

  if(!ifs.good()) {
    ERR << "[PCA::load] Cannot open '"<<filename<<"' to read PCA." << endl;
    return false;
  } else {
    ifs >> dimx >> dimy;
    if(dimx==dimy) dim_=dimx;
    else {
      ERR << "[PCA::load] Strange: PCA not square: dimx="<< dimx << " dimy=" << dimy << endl;
      return false;
    }
    ::std::string tmpstr;
    ifs >> tmpstr;
    mean_=vector<double>(dimy);
    for(int i=0;i<dimy;++i) {
      ifs >> mean_[i];
    }
    covariance_=vector< vector<double> >(dimy,vector<double>(dimx,0.0));
    for(int y=0;y<dimy;++y) {
      ifs >> posy;
      if(posy!=y) {
	ERR << "[PCA::load] Strangeness in reading PCA." << endl;
	return false;
      }
      for(int x=0;x<dimx;++x) {
        ifs >> covariance_[y][x];
      }
    }
    ifs >> posy;
    if(posy!=-1) {
      ERR << "[PCA::load] Reading PCA was strange. Did not end with -1." << endl;
      return false;
    }
  }
  return true;
}

void PCA::dataEnd() {
  for(int i=0;i<dim_;++i) {
    mean_[i]/=counter_;
  }
  
  for(int i=0;i<dim_;++i) {
    for(int j=0;j<dim_;++j) {
      covariance_[i][j]/=counter_;
    }
  }
  
  for(int i=0;i<dim_;++i) {
    for(int j=0;j<dim_;++j) {
      covariance_[i][j]-=mean_[i]*mean_[j];
    }
  }
}

void PCA::calcPCA() {
  
  vector<double> w(dim_);
  vector< vector<double> > v(dim_, vector<double>(dim_,0));

#ifdef HAVE_LAPACK
  svd(covariance_,w,v);
  //  NR::svdcmp(covariance_, w, v);
#else
#warning "You are compiling PCA without SVD. This won't work"
  ERR << "PCA running without SVD. " << endl
      << "some things after this point are definitively broken." << endl;
#endif
      
  DBG(50) << "Eigenvalues: " ;
  for(unsigned int i=0;i<w.size();++i) {
    BLINK(50) << w[i] << " ";
  }
  BLINK(50) << endl;

  DBGI(50,{
      DBG(50) << "transform before sorting:" << endl;
      for(int i=0;i<dim_;++i) {
        for(int j=0;j<dim_;++j) {
          BLINK(50) << " " << covariance_[i][j] ;
        }
        BLINK(50) << endl;
      }
    });

  vector< pair<double, vector<double> > > toSort;
  for(unsigned int i=0;i<w.size();++i) {
    vector<double> tmp(dim_);
    for(int j=0;j<dim_;++j) tmp[j]=covariance_[j][i];
    toSort.push_back(pair<double,vector<double> >(w[i],tmp));
  }
  
  sort(toSort.rbegin(),toSort.rend());
  
  for(unsigned int i=0;i<toSort.size();++i) {
    covariance_[i]=toSort[i].second;
    w[i]=toSort[i].first;
  }
  DBG(50) << "Eigenvalues: " ;
  for(unsigned int i=0;i<w.size();++i) {
    BLINK(50) << w[i] << " ";
  }
  BLINK(50) << endl;
  
  DBG(50) << "transform after sorting:" << endl;
  for(int i=0;i<dim_;++i) {
    for(int j=0;j<dim_;++j) {
      BLINK(50) << " " << covariance_[i][j] ;
    }
    BLINK(50)  << endl;
  }

  eigenvalues_=w;
  
}

const int PCA::dim() const {
  return dim_;
}

const vector<double> & PCA::mean() const {
  return mean_;
}

const vector< vector <double> > & PCA::covariance() const {
  return covariance_;
}

const int PCA::counter() const {
  return counter_;
}

const vector<double> PCA::transform(const vector<double> &in, int dim) const {
  if(dim==0) dim=in.size();
  vector<double> result(dim,0.0);
  DBG(20) << "Transforming to dim " << dim << endl;
  
  for(int i=0;i<dim;++i) {
    for(unsigned int j=0;j<in.size();++j) {
      result[i]+=covariance_[i][j]*(in[j]-mean_[j]);
    }
  }
  return result;
}



/** as PCA-matrices are rotation matrices, for this calculation it is
 *  not necessary to explicitly invert them, but instead it holds:
 *
 *   M^{-1}=transpose(M), and thus M×transpose(M)=Id
 */
const vector<double> PCA::backTransform(vector<double> in) const {
  if(int(in.size()) < dim_) { in.resize(dim_,0);} // zeropadding
  vector<double> result(dim_,0.0);
  DBG(20) << "Backtransforming" << endl;
  
  for(int i=0;i<dim_;++i) {
    for(int j=0;j<dim_;++j) {
      result[i]+=covariance_[j][i]*in[j];
    }
    result[i]+=mean_[i];
  }
  return result;
}
