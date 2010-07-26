#include "dist_weightedl1.hpp"
#include <iostream>
#include "diag.hpp"


using namespace std;

double WeightedL1Distance::distance(const BaseFeature* queryFeature, const BaseFeature* databaseFeature) {
  double result=0.0;
  const VectorFeature* db=dynamic_cast<const VectorFeature*>(databaseFeature);
  const VectorFeature* query=dynamic_cast<const VectorFeature*>(queryFeature);
  if (db && query) {
    double tmp;
    for (uint i=0; i<db->size(); ++i) {
      tmp=(db->operator[](i))-(query->operator[](i));
      result+=weight_[i]*fabs(tmp);
    }
    return result;

  } else {
    ERR << "Features not comparable: need vector features" << ::std::endl;
    return -1.0;
  }
}



void WeightedL1Distance::tune(const std::vector<const BaseFeature*>& posFeat, const std::vector<const BaseFeature*>& negFeat) {

  int dim=dynamic_cast<const VectorFeature*>(posFeat[0])->size();
  weight_=vector<double>(dim,1.0);
  if(negFeat.size()<2 or posFeat.size()<2) {
    DBG(10) << "Not enough data to train weights. Returning" << endl;
    DBG(50) << "weights:";
    return;
  }
  
  int i,j,k,it;
  int jdn=0,jnn=0;  
  float index,dn,dd,fv,dif,dnk,ddk,dis;
  float expon,sigmoid;
  float reg;
  int p,n=0;

  int samples=posFeat.size()+negFeat.size();
  float *datos[samples];
  int datclas[samples];
  int sel[samples];
 
  DBG(50) << VAR(dim) << " " << VAR(samples) << endl;
  for(int d=0;d<dim;++d) {
    BLINK(50) << " " <<weight_[d];
  }
  BLINK(50) << endl;

  for(uint no=0;no<posFeat.size();++n,++no) {
    const VectorFeature* vec=dynamic_cast<const VectorFeature*>(posFeat[no]);
    if(!vec) {ERR << "Cannot tune the distance. Need Vector Features" << endl; return;}
    datos[n]=new float[dim];
    
    for(int d=0;d<dim;++d) {
      datos[n][d]=(*vec)[d];
    }
    sel[n]=1;
    datclas[n]=1;
  }

  for(uint no=0;no<negFeat.size();++n,++no) {
    const VectorFeature* vec=dynamic_cast<const VectorFeature*>(negFeat[no]);
    if(!vec) {ERR << "Cannot tune the distance. Need Vector Features" << endl; return;}
    datos[n]=new float[dim];
    for(int d=0;d<dim;++d) {
      datos[n][d]=(*vec)[d];
    }
    sel[n]=-1;
    datclas[n]=0;
  }
  
  float sg[dim];
  for(int d=0;d<dim;++d) {
    sg[d]=1.0; weight_[d]=1.0;
  }
  
  p=n=0;
  for(i=0;i<samples;i++) {
    if (sel[i]==1) p++;
  }
    
  for(i=0;i<samples;i++) {    
      if (sel[i]==-1) n++;
  }

  if ((n==0)||(p<2)) return;

     
  it=0;
  while (it<maxIterations_) {
    it++;
    index=0.0;

    for (i=0; i<samples; i++) {
      if (sel[i]==1) {
        dn=FLT_MAX;
        for (j=0; j<samples; j++)
          if ((j!=i)&&(sel[j]==1)) {
            dis=wdisl1(datos[i], datos[j], dim);
            if (dis<dn) {
              dn=dis;
              jnn=j;
            }
          }

        dd=FLT_MAX;
        for (j=0; j<samples; j++)
          if (sel[j]==-1) {
            dis=wdisl1(datos[i], datos[j], dim);
            if (dis<dd) {
              dd=dis;
              jdn=j;
            }
          }

        index+=sigmoide(dn/dd);

        //gradient
        expon=exp(rampa-(rampa*(dn/dd)));
        sigmoid=(rampa*expon)/((1+expon)*(1+expon));
        for (k=0; k<dim; k++) {// Update sigmas
          dnk=weight_[k]*fabs(datos[jnn][k]-datos[i][k]);
          ddk=weight_[k]*fabs(datos[jdn][k]-datos[i][k]);
          if ((dnk!=0.0)&&(ddk!=0.0)) {
            dif=(datos[jnn][k]-datos[i][k])*(datos[jnn][k]-datos[i][k]);
            fv=weight_[k]*dif*dd*ddk;
            dif=(datos[jdn][k]-datos[i][k])*(datos[jdn][k]-datos[i][k]);
            fv-=weight_[k]*dif*dn*dnk;
            fv/=ddk*dnk*dd*dd;
            fv*=sigmoid;
            // regularisation
            reg=(1.0-sg[k]);
            sg[k]-=gradientWeight_*(stepWidth_*fv)/p - regularisationWeight_*reg;
            if (sg[k]<0.0)
              sg[k]=0.00001;
          }
        }

      }
    }

    for (i=0; i<samples; i++) {
      if (sel[i]==-1) {

        dn=FLT_MAX;
        for (j=0; j<samples; j++)
          if (sel[j]==1) {
            dis=wdisl1(datos[i], datos[j], dim);
            if (dis<dn) {
              dn=dis;
              jnn=j;
            }
          }

        dd=FLT_MAX;
        for (j=0; j<samples; j++)
          if ((sel[j]==-1)&&(j!=i)) {
            dis=wdisl1(datos[i], datos[j], dim);
            if (dis<dd) {
              dd=dis;
              jdn=j;
            }
          }

        index+=sigmoide(dn/dd);

        //gradient
        expon=exp(rampa-(rampa*(dn/dd)));
        sigmoid=(rampa*expon)/((1+expon)*(1+expon));
        for (k=0; k<dim; k++) {// Update sigmas
          dnk=weight_[k]*fabs(datos[jnn][k]-datos[i][k]);
          ddk=weight_[k]*fabs(datos[jdn][k]-datos[i][k]);
          if ((dnk!=0.0)&&(ddk!=0.0)) {
            dif=(datos[jnn][k]-datos[i][k])*(datos[jnn][k]-datos[i][k]);
            fv=weight_[k]*dif*dd*ddk;

            dif=(datos[jdn][k]-datos[i][k])*(datos[jdn][k]-datos[i][k]);
            fv-=weight_[k]*dif*dn*dnk;

            fv/=ddk*dnk*dd*dd;
            fv*=sigmoid;
            //regularistion
            reg=(1.0-sg[k]);
            sg[k]+=gradientWeight_*(stepWidth_*fv)/n + regularisationWeight_*reg;
            if (sg[k]<0.0)
              sg[k]=0.00001;
          }
        }
      }
    }
    index/=(p+n);
    DBG(50) << VAR(it) << " " << VAR(index) << "weights:";    
    for(int d=0;d<dim;++d) {
      weight_[d]=sg[d];
      BLINK(50) << " " <<weight_[d];
    }
    BLINK(50) << endl;
    
  }

  DBG(50) << "weights:";
  for(int d=0;d<dim;++d) {
    BLINK(50) << " " <<weight_[d];
  }
  BLINK(50) << endl;

  for(int n=0;n<samples;++n) {delete[] datos[n];}
}

float WeightedL1Distance::wdisl1(float * i, float *j, int dim) {
  float dis=0.0;
  for(int d=0;d<dim;++d) {
    dis+=weight_[d]*fabs(i[d]-j[d]);
  }
  return dis;
}


float WeightedL1Distance::sigmoide(float x){
  return 1/(1+exp(rampa-(rampa*x)));
}
