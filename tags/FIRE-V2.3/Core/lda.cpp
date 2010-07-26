#include "lda.hpp"
#include "Lapack.hh"
#include "diag.hpp"
#include "gzstream.hpp"
#include <limits>

using namespace std;
using namespace Lapack;

LDA::LDA(int dim, int noCls) : dim_(dim), noCls_(noCls), counter_(0),
                               withinScatter_(dim,vector<double>(dim)), 
                               betweenScatter_(dim,vector<double>(dim)),
                               eigenVecs_(dim,vector<double>(dim)),
                               classMeans_(noCls,vector<double>(dim)),
                               eigenVals_(dim),
                               mean_(dim),
                               observationsInClasses_(noCls)
{
  
}

void LDA::putDataFirst(const ::std::vector<double> &in, int clasz) {

  ++observationsInClasses_[clasz];
  ++counter_;
  for(int d=0;d<dim_;++d) {
    classMeans_[clasz][d]+=in[d];
  }
}

void LDA::dataEndFirst() {
  for(int c=0;c<noCls_;++c) {
    for(int d=0;d<dim_;++d) {
      classMeans_[c][d]/=double(observationsInClasses_[c]);
      mean_[d]+=classMeans_[c][d];
    }
  }
  for(int d=0;d<dim_;++d) {
    mean_[d]/=double(noCls_);
  }
}


void LDA::putDataSecond(const ::std::vector<double> &in, int clasz) {

  // within Scatter
  for(int d=0;d<dim_;++d) {
    for(int t=0;t<=d;++t) {
      withinScatter_[d][t]+=(in[d]-classMeans_[clasz][d])*(in[t]-classMeans_[clasz][t]);
    }
  }
  
  // between Scatter
  for(int c=0;c<noCls_;++c) {
    for(int d=0;d<dim_;++d) {
      for(int t=0;t<=d;++t) {
        betweenScatter_[d][t]+=double(observationsInClasses_[c])*(classMeans_[c][d]-mean_[d])*(classMeans_[c][t]-mean_[t]);
      }
    }
  }
}

void LDA::dataEndSecond() {
  for(int d=0;d<dim_;++d) {
    for(int t=0;t<=d;++t) {
      betweenScatter_[d][t]/=double(counter_);
      betweenScatter_[t][d]=betweenScatter_[d][t];
      withinScatter_[d][t]/=double(counter_);
      withinScatter_[t][d]=withinScatter_[d][t];
    }
  }
}


void LDA::calcLDA() {
  char * jobYes="V";
  char * jobNo="N";
  
  int n_dim_integer=dim_;
  double *between_scatter=new double[dim_*dim_];
  double *within_scatter=new double[dim_*dim_];
  double *eigenvals_r=new double[dim_];
  double *eigenvals_i=new  double[dim_];
  double *eigenvals_b=new double[dim_];
  double *eigenvecs=new double[dim_*dim_];
  int lwork=(8+dim_)*dim_;
  double *work=new double[lwork];
  int info;

  for(int d=0;d<dim_;++d) {
    for(int t=0;t<dim_;++t) {
      between_scatter[d+t*dim_]=betweenScatter_[d][t];
      within_scatter[d+t*dim_]=withinScatter_[d][t];
    }
  }

  DBGI(20,{
      cout << "between_scatter " << endl;
      for(int d=0;d<dim_;++d) {
        for(int t=0;t<dim_;++t) {
          cout << between_scatter[d+t*dim_] << " ";
        }
        cout << endl;
      }
      
      cout << endl << "within_scatter" << endl;
      for(int d=0;d<dim_;++d) {
        for(int t=0;t<dim_;++t) {
          cout << within_scatter[d+t*dim_] << " ";
        }
        cout << endl;
      }
    });
  
  DBG(20) << "calling lapack..." << endl;

  dggev(jobYes, jobYes,
        &n_dim_integer,
        between_scatter, &n_dim_integer, 
        within_scatter,  &n_dim_integer, 
        eigenvals_r, eigenvals_i, eigenvals_b,
        eigenvecs, &n_dim_integer,
        eigenvecs, &n_dim_integer,
        work, &lwork, &info);

  DBG(10) << "lapack answered." << endl;

  if(info==0) {
    DBG(10) << "LDA successfully finished" << endl;
  } else {
    ERR << "LDA (sggev) returned " << info << endl;
    exit(20);
  }

  DBGI(20,{
      cout << "eigenvals_r"<< endl;
      for(int d=0;d<dim_;++d) {
        cout << eigenvals_r[d] << " ";
      }
      cout << endl;

      cout<< endl << "eigenvals_i"<< endl;
      for(int d=0;d<dim_;++d) {
        cout << eigenvals_i[d] << " ";
      }
      cout << endl;

      cout<< endl << "eigenvals_b"<< endl;
      for(int d=0;d<dim_;++d) {
        cout << eigenvals_b[d] << " ";
      }
      cout << endl;
    });
  
  // copy eigenvalues
  DBG(20) << "eigenvalues " ;
  for(int d=0;d<dim_;++d) {
    if(eigenvals_b[d]<numeric_limits<double>::epsilon()) {
      eigenVals_[d]=-numeric_limits<double>::max();
    } else {
      eigenVals_[d]=eigenvals_r[d]/eigenvals_b[d];
    }
    BLINK(20) << eigenVals_[d] << " ";
  }
  

  vector< pair<double, vector<double> > > toSort;
  for(int i=0;i<dim_;++i) {
    vector<double> tmp(dim_);
    for(int d=0;d<dim_;++d) {
      tmp[d]=eigenvecs[d+i*dim_];
    }
    toSort.push_back( pair<double, vector<double> >(eigenVals_[i],tmp));
  }
  sort(toSort.rbegin(), toSort.rend());
    
  BLINK(20)<< endl;
  // copy eigenvectors
  for(int d=0;d<dim_;++d) {
    eigenVals_[d]=toSort[d].first;
    for(int t=0;t<dim_;++t) {
      eigenVecs_[d][t]=toSort[d].second[t];
    }
  }

  
  delete[] between_scatter;
  delete[] within_scatter;
  delete[] eigenvals_r;
  delete[] eigenvals_i;
  delete[] eigenvals_b;
  delete[] eigenvecs;
  delete[] work;
}


const vector<double> LDA::transform(const vector<double> &in, int dim) const {
  if(dim==0) dim=in.size();
  vector<double> result(dim,0.0);
  DBG(20) << "Transforming to dim " << dim << endl;
  
  for(int i=0;i<dim;++i) {
    for(unsigned int j=0;j<in.size();++j) {
      result[i]+=eigenVecs_[i][j]*(in[j]-mean_[j]);
    }
  }
  return result;
}

/** as LDA-matrices are rotation matrices, for this calculation it is
 *  not necessary to explicitly invert them, but instead it holds:
 *
 *   M^{-1}=transpose(M), and thus M×transpose(M)=Id
 */
const vector<double> LDA::backTransform(vector<double> in) const {
  if(int(in.size()) < dim_) { in.resize(dim_,0);} // zeropadding
  vector<double> result(dim_,0.0);
  DBG(20) << "Backtransforming" << endl;
  
  for(int i=0;i<dim_;++i) {
    for(int j=0;j<dim_;++j) {
      result[i]+=eigenVecs_[j][i]*in[j];
    }
    result[i]+=mean_[i];
  }
  return result;
}



void LDA::save(const string &filename) const{
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
        ofs << " " << eigenVecs_[y][x];
      }
      ofs << endl;
    }
    ofs << -1 << endl;
  }
}

void LDA::load(const string &filename) {
  igzstream ifs; ifs.open(filename.c_str());
  int dimy, dimx, posy;

  if(!ifs.good()) {
    ERR << "Cannot open '"<<filename<<"' to read PCA." << endl;
  } else {
    ifs >> dimx >> dimy;
    if(dimx==dimy) dim_=dimx;
    else ERR << "Strange: LDA not square: dimx="<< dimx << " dimy=" << dimy << endl;
    ::std::string tmpstr;
    ifs >> tmpstr;
    mean_=vector<double>(dimy);
    for(int i=0;i<dimy;++i) {
      ifs >> mean_[i];
    }
    eigenVecs_=vector< vector<double> >(dimy,vector<double>(dimx,0.0));
    for(int y=0;y<dimy;++y) {
      ifs >> posy;
      if(posy!=y) {ERR << " Strangeness in reading." << endl;}
      for(int x=0;x<dimx;++x) {
        ifs >> eigenVecs_[y][x];
      }
    }
    ifs >> posy;
    if(posy!=-1) {ERR << "Reading was strange. Did not end with -1." << endl;}
  }
}
