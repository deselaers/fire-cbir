#include "querycombiner.hpp"
#include "stringparser.hpp"
#include "supportvectormachine.hpp"
#include "vectorfeature.hpp"
#include <sstream>
#include <vector>
using namespace std;

ScoreSumQueryCombiner::ScoreSumQueryCombiner(Retriever& retriever) :
  posWeight_(1.0), negWeight_(0.8333), retriever_(retriever) {
}

ScoreSumQueryCombiner::~ScoreSumQueryCombiner() {
}

void ScoreSumQueryCombiner::query(const vector<ImageContainer*> posQueries, const vector<ImageContainer*> negQueries, vector<ResultPair>& results) {

  double posWeight=posWeight_/posQueries.size();
  double negWeight=negWeight_/negQueries.size();
  uint N=retriever_.database().size();
  retriever_.imageComparator().tuneDistances(posQueries,negQueries);
  
  vector<double> activeScores(N, 0.0);
  DBG(10) << "Querying ";
  for (long q=0; q<long(posQueries.size()); ++q) {
    BLINK(10) << "+";
    DBG(15) << "Positive query: " << posQueries[q]->basename() << endl;
    if (retriever_.imageComparator().size() != posQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than positive query " << q << " (" << posQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(posQueries[q]);
    retriever_.getScores(posQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      results[i].first+=posWeight*activeScores[i];
    }
  }

  for (long q=0; q<long(negQueries.size()); ++q) {
    BLINK(10) << "-";
    DBG(15) << "Negative query: " << negQueries[q]->basename() << endl;

    if (retriever_.imageComparator().size() != negQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than negative query " << q << " (" << negQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(negQueries[q]);
    retriever_.getScores(negQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      results[i].first+=negWeight*(1.0-activeScores[i]);
    }

  }
  BLINK(10) << endl;
}

void ScoreSumQueryCombiner::setParameters(const std::string & parameters) {
  posWeight_=getDoubleAfter(parameters, "POS=", 1.0);
  negWeight_=getDoubleAfter(parameters, "NEG=", 0.8333);
}

NNQuotientQueryCombiner::NNQuotientQueryCombiner(Retriever& retriever) :
  retriever_(retriever) {
}

NNQuotientQueryCombiner::~NNQuotientQueryCombiner() {

}

void NNQuotientQueryCombiner::query(const vector<ImageContainer*> posQueries, const vector<ImageContainer*> negQueries, vector<ResultPair>& results) {

  uint N=retriever_.database().size();
  vector<double> activeScores(N, 0.0), bestPos(N, 0.0), bestNeg(N, 0.0);
  
  retriever_.imageComparator().tuneDistances(posQueries,negQueries);
  
  DBG(10) << "Querying ";
  // first get scores for all positive queries
  for (long q=0; q<long(posQueries.size()); ++q) {
    DBG(15) << "Positive query: " << posQueries[q]->basename() << endl;
    if (retriever_.imageComparator().size() != posQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than positive query " << q << " (" << posQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(posQueries[q]);
    retriever_.getScores(posQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      if (bestPos[i]<activeScores[i]) {
        bestPos[i]=activeScores[i];
      }
    }
  }

  if (negQueries.size()==0) {
    bestNeg=vector<double>(N, 1.0);
  }
  // then get scores fr all negative queries
  for (long q=0; q<long(negQueries.size()); ++q) {
    BLINK(10) << "-";
    DBG(15) << "Negative query: " << negQueries[q]->basename() << endl;
    if (retriever_.imageComparator().size() != negQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than negative query " << q << " (" << negQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(negQueries[q]);
    retriever_.getScores(negQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      if (bestNeg[i]<activeScores[i]) {
        bestNeg[i]=activeScores[i];
      }
    }
  }
  BLINK(10) << endl;

  // now combine these scores
  for (uint n=0; n<N; ++n) {
    results[n].first=bestPos[n]/bestNeg[n];
  }

}

SVMQueryCombiner::SVMQueryCombiner(Retriever& retriever) :
  retriever_(retriever) {
}

SVMQueryCombiner::~SVMQueryCombiner() {
}

void SVMQueryCombiner::query(const vector<ImageContainer*> posQueries, const vector<ImageContainer*> negQueries, vector<ResultPair>& results) {
  vector<DoubleVector> trainingExamples;
  vector<int> classes;

  if (negQueries.size()==0) {
    DBG(10) << "Only one example: NN" << endl;
    ScoreSumQueryCombiner adc(retriever_);
    adc.query(posQueries, negQueries, results);
  } else {
    DBG(10) << "Training new SVM" << endl;
    trainingExamples.reserve(posQueries.size()+negQueries.size());
    classes.reserve(posQueries.size()+negQueries.size());

    for (uint q=0; q<posQueries.size(); ++q) {
      trainingExamples.push_back(posQueries[q]->asVector());
      classes.push_back(1);
    }
    for (uint q=0; q<negQueries.size(); ++q) {
      trainingExamples.push_back(negQueries[q]->asVector());
      classes.push_back(0);
    }
    SupportVectorMachine svm(kerneltype_, degree_, gamma_, coef0_, cost_);
    DBG(10) << "training SVM with " << trainingExamples.size() << " vectors of size " << trainingExamples[0].size() <<endl;
    svm.train(trainingExamples, classes);
    uint N=retriever_.database().size();
    DoubleVector scores;
    for (uint n=0; n<N; ++n) {
      svm.classify(retriever_.database()[n]->asVector(), scores);
      results[n].first=scores[0];
    }
  }
}

void SVMQueryCombiner::setParameters(const std::string & parameters) {
  kerneltype_=getIntAfter(parameters, "KERNEL=", 2);
  degree_=getIntAfter(parameters, "DEGREE=", 3);
  gamma_=getDoubleAfter(parameters, "GAMMA=", 0.01);
  coef0_=getDoubleAfter(parameters, "COEF0=", 0.0);
  cost_=getDoubleAfter(parameters, "COST=", 4.0);
}

SumQuotientQueryCombiner::SumQuotientQueryCombiner(Retriever& retriever) :
  retriever_(retriever) {
}

SumQuotientQueryCombiner::~SumQuotientQueryCombiner() {
}

void SumQuotientQueryCombiner::query(const vector<ImageContainer*> posQueries, const vector<ImageContainer*> negQueries, vector<ResultPair>& results) {

  uint N=retriever_.database().size();

  vector<double> activeScores(N, 0.0), posScores(N, 0.0), negScores(N, 0.0);
  
  retriever_.imageComparator().tuneDistances(posQueries,negQueries);
  
  DBG(10) << "Querying ";
  for (long q=0; q<long(posQueries.size()); ++q) {
    DBG(15) << "Positive query: " << posQueries[q]->basename() << endl;
    BLINK(10) << "+";

    if (retriever_.imageComparator().size() != posQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than positive query " << q << " (" << posQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(posQueries[q]);
    retriever_.getScores(posQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      posScores[i]+=activeScores[i];
    }
  }

  if (negQueries.size()==0) {
    negScores=vector<double>(N, 1.0);
  }
  for (long q=0; q<long(negQueries.size()); ++q) {
    DBG(15) << "Negative query: " << negQueries[q]->basename() << endl;
    BLINK(10) << "-";
    if (retriever_.imageComparator().size() != negQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than negative query " << q << " (" << negQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(negQueries[q]);
    retriever_.getScores(negQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      negScores[i]+=activeScores[i];
    }
  }
  BLINK(10) << endl;

  for (uint n=0; n<N; ++n) {
    results[n].first=posScores[n]/negScores[n];
  }
}

RelevanceScoreQueryCombiner::RelevanceScoreQueryCombiner(Retriever& retriever) :
  retriever_(retriever) {
}

RelevanceScoreQueryCombiner::~RelevanceScoreQueryCombiner() {

}

void RelevanceScoreQueryCombiner::query(const vector<ImageContainer*> posQueries, const vector<ImageContainer*> negQueries, vector<ResultPair>& results) {

  uint N=retriever_.database().size();
  vector<double> activeScores(N, 0.0), bestPos(N, std::numeric_limits<double>::epsilon()), bestNeg(N, std::numeric_limits<double>::epsilon());

  retriever_.imageComparator().tuneDistances(posQueries,negQueries);
  
  // first get scores for all positive queries
  DBG(10) << "Querying ";
  for (long q=0; q<long(posQueries.size()); ++q) {
    BLINK(10) << "+";
    DBG(15) << "Positive query: " << posQueries[q]->basename() << endl;
    if (retriever_.imageComparator().size() != posQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than positive query " << q << " (" << posQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(posQueries[q]);
    retriever_.getScores(posQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      if (bestPos[i]<activeScores[i]) {
        bestPos[i]=activeScores[i];
      }
    }
  }

  if (negQueries.size()==0) {
    bestNeg=vector<double>(N, 1e-200);
  }
  // then get scores fr all negative queries
  for (long q=0; q<long(negQueries.size()); ++q) {
    DBG(15) << "Negative query: " << negQueries[q]->basename() << endl;
    BLINK(10) << "-";
    if (retriever_.imageComparator().size() != negQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than negative query " << q << " (" << negQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(negQueries[q]);
    retriever_.getScores(negQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      if (bestNeg[i]<activeScores[i]) {
        bestNeg[i]=activeScores[i];
      }
    }
  }
  BLINK(10) << endl;
  // now combine these scores
  for (uint n=0; n<N; ++n) {
    results[n].first=1.0/(1.0+log(bestPos[n])/log(bestNeg[n]));
  }

}

DistSumQuotientQueryCombiner::DistSumQuotientQueryCombiner(Retriever& retriever) :
  retriever_(retriever) {
}

void DistSumQuotientQueryCombiner::setParameters(const std::string &parameters) {
  string modestring=getStringAfter(parameters,"MODE=","NEGDENOM");
  if(modestring=="NEGDENOM") {
    DBG(10) << "NEGDENOM" << endl;
    mode_=negDenom;
  } else if(modestring=="ALLDENOM") {
    DBG(10) << "ALLDENOM" << endl;
    mode_=allDenom;
  } else {
    ERR << "Unknown setting. Going default: NEGDEMON";
    mode_=negDenom;
  }
}


DistSumQuotientQueryCombiner::~DistSumQuotientQueryCombiner() {
}

void DistSumQuotientQueryCombiner::query(const vector<ImageContainer*> posQueries, const vector<ImageContainer*> negQueries, vector<ResultPair>& results) {
  uint N=retriever_.database().size();

  vector<double> activeScores(N, 0.0), posScores(N, 0.0), negScores(N, 0.0);

  retriever_.imageComparator().tuneDistances(posQueries,negQueries);
  
  DBG(10) << "Querying ";
  for (long q=0; q<long(posQueries.size()); ++q) {
    DBG(15) << "Positive query: " << posQueries[q]->basename() << endl;
    BLINK(10) << "+";
    if (retriever_.imageComparator().size() != posQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than positive query " << q << " (" << posQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(posQueries[q]);
    retriever_.getScores(posQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      posScores[i]+=1.0/(log(activeScores[i])+0.0001);
    }
  }

  if (negQueries.size()==0) {
    negScores=vector<double>(N, 1.0);
  }
  for (long q=0; q<long(negQueries.size()); ++q) {
    DBG(15) << "Negative query: " << negQueries[q]->basename() << endl;
    BLINK(10) << "-";
    if (retriever_.imageComparator().size() != negQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than negative query " << q << " (" << negQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(negQueries[q]);
    retriever_.getScores(negQueries[q], activeScores);
    retriever_.imageComparator().stop();
    for (uint i=0; i<N; ++i) {
      negScores[i]+=1.0/(log(activeScores[i])+0.0001);
    }
  }
  BLINK(10) << endl;
  if(mode_==negDenom) {
    for (uint n=0; n<N; ++n) {
      results[n].first=1.0/(1.0+posScores[n]/negScores[n]);
    }
  } else {
    for (uint n=0; n<N; ++n) {
      results[n].first=(posScores[n])/(negScores[n]+posScores[n]);
    }
  }
}


WeightedDistanceQueryCombiner::WeightedDistanceQueryCombiner(Retriever& retriever) :
  retriever_(retriever), rampa(10.0), FLT_MAX(100000000.0), FLT_MIN(-100000000.0), mode_(WDRelevanceScore) {
}

WeightedDistanceQueryCombiner::~WeightedDistanceQueryCombiner() {
}

void WeightedDistanceQueryCombiner::query(const vector<ImageContainer*> posQueries, const vector<ImageContainer*> negQueries, vector<ResultPair>& results) {
  
  if (negQueries.size()==0) {
    DBG(10) << "Only positive examples: defaulting to ScoreSumQueryCombiner" << endl;
    ScoreSumQueryCombiner adc(retriever_);
    adc.query(posQueries, negQueries, results);
  } else {
    DBG(10) << "Training new WeightedDistance" << endl;
    
    samples=posQueries.size()+negQueries.size();
    datos=new float*[samples];
    datclas=new int[samples];
    sel=new int[samples];
    dim=posQueries[0]->asVector().size();
    sigmas=new float[dim];
    classes=2;
    
    // init sigmas:
    for(int i=0;i<dim;++i) sigmas[i]=1.0;
    
    // copy query data into appropriate data structures
    DoubleVector tmpVec;
    uint pos=0;
    DBG(10) << "Querying ";
    for (uint q=0; q<posQueries.size(); ++q) {
      BLINK(10) << "+";
      tmpVec=posQueries[q]->asVector();
      datos[pos]=new float[dim];
      for(int d=0;d<dim;++d) {
        datos[pos][d]=tmpVec[d];
      }
      datclas[pos]=1;
      sel[pos]=1;
      ++pos;
    }
    
    for (uint q=0; q<negQueries.size(); ++q) {
      BLINK(10) << "-";
      tmpVec=negQueries[q]->asVector();
      datos[pos]=new float[dim];
      for(int d=0;d<dim;++d) {
        datos[pos][d]=tmpVec[d];
      }
      datclas[pos]=0;
      sel[pos]=-1;
      ++pos;
    }
    BLINK(10) << endl;
    calcWeights();
    DBG(10) << VAR(printSigmas()) << endl;

    // now get retrieval result: calculate a score for each image in the database
    uint N=retriever_.database().size();
    switch(mode_) {
    case WDRelevanceScore:
      for (uint n=0; n<N; ++n) {
        double minRel=std::numeric_limits<double>::max(),  minNRel=std::numeric_limits<double>::max();
        for(int i=0;i<samples;++i) {
          double dis=WeightedL1Distance(datos[i],retriever_.database()[n]->asVector());
          if(sel[i]==1) {
            if(minRel>dis) {minRel=dis;}
          } else if(sel[i]==-1) {
            if(minNRel>dis) {minNRel=dis;}
          }
          results[n].first=minNRel/minRel;
        } 
      }
      break;
    case WDDistSumQuotient:
      for (uint n=0; n<N; ++n) {
        double sumRel=0.0, sumNRel=0.0;
        for(int i=0;i<samples;++i) {
          double dis=WeightedL1Distance(datos[i],retriever_.database()[n]->asVector());
          if(sel[i]==1) {sumRel+=dis;}
          else if(sel[i]==-1) {sumNRel+=dis;}
          else ERR << "Something wrong" << endl;
        }
        results[n].first=sumNRel/sumRel;
      } 
      break;
    case WDScoreSum: {
      // this is old fire's relevance feedback, but with weighted L1 distance
      vector<vector<double> > distances(N,vector<double>(samples));
      vector<double> distsums(samples,0);
      for(uint n=0;n<N;++n) {
        for(int i=0;i<samples;++i) {
          distances[n][i]=WeightedL1Distance(datos[i],retriever_.database()[n]->asVector());
          distsums[i]+=distances[n][i];     
        }
      }
      
      for(uint n=0;n<N;++n) {
        double sumRel=0.0, sumNRel=0.0; int rel=0, nrel=0;
        for(int i=0;i<samples;++i) {
          distances[n][i]=exp(-distances[n][i]/distsums[i]);
          if(sel[i]==1) { sumRel+=distances[n][i]; rel+=1;}
          else if(sel[i]==-1) {sumNRel+=1-distances[n][i]; nrel+=1;}
          else ERR << "Something wrong" <<endl;
        }
        results[n].first=sumRel/rel + 0.8*sumNRel/nrel;
      }
      break;}
    default:
      ERR<< "Unknown mode" << endl;
      break;
    }
   
    // clean up
    for(int n=0;n<samples;++n) {
      delete[] datos[n];
    }
    delete[] datos;
    delete[] datclas;
    delete[] sel;
    delete[] sigmas;
  }
}

void WeightedDistanceQueryCombiner::copiasig(float *sg, const float *sigmas) {
  int i;
  for (i=0; i<dim; i++)
    sg[i]=sigmas[i];
}
  
void WeightedDistanceQueryCombiner::calcWeights() {
  int i,j,k,it;
  int jdn=0,jnn=0;  
  float index,dn,dd,fv,dif,dnk,ddk,dis;
  float stepWidth_=.1,*sg,expon,sigmoid;
  float reg;
  int p,n;

  p=n=0;
  for(i=0;i<samples;i++) {
    if (sel[i]==1) p++;
  }
    
  for(i=0;i<samples;i++) {    
      if (sel[i]==-1) n++;
  }

  if ((n==0)||(p<2)) return;

  sg=(float *)malloc(dim*sizeof(float));
  copiasig(sg,sigmas);
   
  it=0;
  while (it<maxIterations_) {
    it++;
    index=0.0;

    for (i=0; i<samples; i++) {
      if (sel[i]==1) {
        dn=FLT_MAX;
        for (j=0; j<samples; j++)
          if ((j!=i)&&(sel[j]==1)) {
            dis=wdisl1(i, j);
            if (dis<dn) {
              dn=dis;
              jnn=j;
            }
          }

        dd=FLT_MAX;
        for (j=0; j<samples; j++)
          if (sel[j]==-1) {
            dis=wdisl1(i, j);
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
          dnk=sigmas[k]*fabs(datos[jnn][k]-datos[i][k]);
          ddk=sigmas[k]*fabs(datos[jdn][k]-datos[i][k]);
          if ((dnk!=0.0)&&(ddk!=0.0)) {
            dif=(datos[jnn][k]-datos[i][k])*(datos[jnn][k]-datos[i][k]);
            fv=sigmas[k]*dif*dd*ddk;

            dif=(datos[jdn][k]-datos[i][k])*(datos[jdn][k]-datos[i][k]);
            fv-=sigmas[k]*dif*dn*dnk;

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
            dis=wdisl1(i, j);
            if (dis<dn) {
              dn=dis;
              jnn=j;
            }
          }

        dd=FLT_MAX;
        for (j=0; j<samples; j++)
          if ((sel[j]==-1)&&(j!=i)) {
            dis=wdisl1(i, j);
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
          dnk=sigmas[k]*fabs(datos[jnn][k]-datos[i][k]);
          ddk=sigmas[k]*fabs(datos[jdn][k]-datos[i][k]);
          if ((dnk!=0.0)&&(ddk!=0.0)) {
            dif=(datos[jnn][k]-datos[i][k])*(datos[jnn][k]-datos[i][k]);
            fv=sigmas[k]*dif*dd*ddk;

            dif=(datos[jdn][k]-datos[i][k])*(datos[jdn][k]-datos[i][k]);
            fv-=sigmas[k]*dif*dn*dnk;

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

    copiasig(sigmas, sg);
    //}
    //fprintf(stderr,"%d %d : %f\n",n,p,index);
    //if (it==1) fprintf(stderr,"%d %d: index %f",p,n,index);
    //if (it==maxit) fprintf(stderr," --> %f\n",index);
  }
  
  free(sg);
}
  
double WeightedDistanceQueryCombiner::WeightedL1Distance(float *v1, const vector<double>& v2) {
  double result=0.0;
  for(int d=0;d<dim;++d) {
    result+=sigmas[d]*fabs(v1[d]-v2[d]);
  }
  return result;
}

void WeightedDistanceQueryCombiner::setParameters(const std::string & parameters) {
    stepWidth_=getDoubleAfter(parameters,"STEPWIDTH=",0.1);
    maxIterations_=getIntAfter(parameters,"MAXITER=",100);
    regularisationWeight_=getDoubleAfter(parameters,"REGWEIGHT=",0.0);
    gradientWeight_=1.0-regularisationWeight_;
  
    string modestring=getStringAfter(parameters,"MODE=","WDRelevanceScore");
    if(modestring=="WDRelevanceScore") {
      DBG(10) << "mode=WDRelevanceScore" << endl;
      mode_=WDRelevanceScore;
    } else if(modestring=="WDDistSumQuotient") {
      DBG(10) << "mode=WDDistSumQuotient"<<endl;
      mode_=WDDistSumQuotient;
    } else if (modestring=="WDScoreSum") {
      DBG(10) << "mode==WDScoreSum" <<endl;
      mode_=WDScoreSum;      
    } else {  
      ERR << "Unknown mode. Proceedings with default: WDRelevanceScore" << endl;
      mode_=WDRelevanceScore;
    }
    
    DBG(10) << VAR(stepWidth_) << " " << VAR(maxIterations_) << " " << VAR(regularisationWeight_) << " " << VAR(mode_) << endl;
    
}

float WeightedDistanceQueryCombiner::wdisl1(long int pos1, long int pos2) {
  float dis;
  int i;
  dis=0;
  for (i=0; i<dim; i++) 
    dis+=sigmas[i]*fabs(datos[pos1][i]-datos[pos2][i]);
  return dis;
}

float WeightedDistanceQueryCombiner::sigmoide(float x){
  return 1/(1+exp(rampa-(rampa*x)));
}

string WeightedDistanceQueryCombiner::printSigmas() const {
  ostringstream os; 
  for(int i=0;i<dim;++i) {
    os << sigmas[i] << " ";
  }
  return os.str();
}


ClassDependentWeightedDistanceQueryCombiner::ClassDependentWeightedDistanceQueryCombiner(Retriever& retriever) :
  retriever_(retriever), rampa(10.0), FLT_MAX(100000000.0), FLT_MIN(-100000000.0), mode_(WDRelevanceScore) {
}

ClassDependentWeightedDistanceQueryCombiner::~ClassDependentWeightedDistanceQueryCombiner() {
}

void ClassDependentWeightedDistanceQueryCombiner::query(const vector<ImageContainer*> posQueries, const vector<ImageContainer*> negQueries, vector<ResultPair>& results) {
  
  if (negQueries.size()==0) {
    DBG(10) << "Only positive examples: defaulting to ScoreSumQueryCombiner" << endl;
    ScoreSumQueryCombiner adc(retriever_);
    adc.query(posQueries, negQueries, results);
  } else {
    DBG(10) << "Training new WeightedDistance" << endl;
    
    samples=posQueries.size()+negQueries.size();
    datos=new float*[samples];
    datclas=new int[samples];
    sel=new int[samples];
    dim=posQueries[0]->asVector().size();
    sigmas=new float*[2];
    sigmas[0]=new float[dim]; sigmas[1]=new float[dim];
    classes=2;
    
    // init sigmas:
    for(int j=0;j<2;++j) {for(int i=0;i<dim;++i){ sigmas[j][i]=1.0;}}
    
    // copy query data into appropriate data structures
    DoubleVector tmpVec;
    uint pos=0;
    DBG(10) << "Querying ";
    for (uint q=0; q<posQueries.size(); ++q) {
      BLINK(10) << "+";
      tmpVec=posQueries[q]->asVector();
      datos[pos]=new float[dim];
      for(int d=0;d<dim;++d) {
        datos[pos][d]=tmpVec[d];
      }
      datclas[pos]=1;
      sel[pos]=1;
      ++pos;
    }
    
    for (uint q=0; q<negQueries.size(); ++q) {
      BLINK(10) << "-";
      tmpVec=negQueries[q]->asVector();
      datos[pos]=new float[dim];
      for(int d=0;d<dim;++d) {
        datos[pos][d]=tmpVec[d];
      }
      datclas[pos]=0;
      sel[pos]=-1;
      ++pos;
    }
    BLINK(10) << endl;

    calcWeights();
    DBG(10) << VAR(printSigmas()) << endl;

    // now get retrieval result: calculate a score for each image in the database
    uint N=retriever_.database().size();
    switch(mode_) {
    case WDRelevanceScore:
      for (uint n=0; n<N; ++n) {
        double minRel=std::numeric_limits<double>::max(),  minNRel=std::numeric_limits<double>::max();
        for(int i=0;i<samples;++i) {
          double dis=WeightedL1Distance(datos[i],retriever_.database()[n]->asVector(),datclas[i]);
          if(sel[i]==1) {
            if(minRel>dis) {minRel=dis;}
          } else if(sel[i]==-1) {
            if(minNRel>dis) {minNRel=dis;}
          }
          results[n].first=minNRel/minRel;
        } 
      }
      break;
    case WDDistSumQuotient:
      for (uint n=0; n<N; ++n) {
        double sumRel=0.0, sumNRel=0.0;
        for(int i=0;i<samples;++i) {
          double dis=WeightedL1Distance(datos[i],retriever_.database()[n]->asVector(),datclas[i]);
          if(sel[i]==1) {sumRel+=dis;}
          else if(sel[i]==-1) {sumNRel+=dis;}
        }
        results[n].first=sumNRel/sumRel;
      } 
      break;
    default:
      ERR<< "Unknown mode" << endl;
      break;
    }
   
    // clean up
    for(int n=0;n<samples;++n) {
      delete[] datos[n];
    }
    delete[] datos;
    delete[] datclas;
    delete[] sel;
    delete[] sigmas[0]; delete[] sigmas[1]; delete[] sigmas;
  }
}

void ClassDependentWeightedDistanceQueryCombiner::copiasig(float **sg, float **sigmas) {
  for (int j=0; j<2; ++j) {
    for (int i=0; i<dim; i++) {
      sg[j][i]=sigmas[j][i];
    }
  }
} 

void ClassDependentWeightedDistanceQueryCombiner::calcWeights() {
  int i,j,k,it,maxit=100;
  int jdn=0,jnn=0;
  float index,dn,dd,fv,dif,dnk,ddk,dis;
  float mu=.1,*sg[2],expon,sigmoid;
  int p,n;

  p=n=0;
  for(i=0;i<samples;i++)
  if (sel[i]==1) p++;

  for(i=0;i<samples;i++)
  if (sel[i]==-1) n++;

  if ((n==0)||(p<2)) return;

  sg[0]=(float *)malloc(dim*sizeof(float));
  sg[1]=(float *)malloc(dim*sizeof(float));

  copiasig(sg,sigmas);

  it=0;
  while (it<maxit) {
    it++;
    index=0.0;

    for(i=0;i<samples;i++)
    if (sel[i]==1) {

      dn=FLT_MAX;
      for(j=0;j<samples;j++)
      if ((j!=i)&&(sel[j]==1)) {
        dis=wdisl1(i,j,sel[j]);
        if (dis<dn) {
          dn=dis;
          jnn=j;
        }
      }

      dd=FLT_MAX;
      for(j=0;j<samples;j++)
      if (sel[j]==-1) {
        dis=wdisl1(i,j,sel[j]);
        if (dis<dd) {
          dd=dis;
          jdn=j;
        }
      }

      index+=sigmoide(dn/dd);
      expon=exp(rampa-(rampa*(dn/dd)));
      sigmoid=(rampa*expon)/((1+expon)*(1+expon));
      for(k=0;k<dim;k++) {// Update sigmas
        dnk=sigmas[0][k]*fabs(datos[jnn][k]-datos[i][k]);
        ddk=sigmas[1][k]*fabs(datos[jdn][k]-datos[i][k]);
        if ((dnk!=0.0)&&(ddk!=0.0)) {
          // Relevant sigmas
          dif=(datos[jnn][k]-datos[i][k])*(datos[jnn][k]-datos[i][k]);
          fv=sigmas[0][k]*dif/(dnk*dd);
          fv*=sigmoid;
          sg[0][k]-=(mu*fv)/p;
          if (sg[0][k]<0.0) sg[0][k]=0.00001;

          // Non-Relevant sigmas
          dif=(datos[jdn][k]-datos[i][k])*(datos[jdn][k]-datos[i][k]);
          fv=sigmas[1][k]*dif*dn/(ddk*dd*dd);
          fv*=sigmoid;
          sg[1][k]+=(mu*fv)/n;
          if (sg[1][k]<0.0) sg[1][k]=0.00001;
        }
      }
    }

    for(i=0;i<samples;i++)
    if (sel[i]==-1) {
      dn=FLT_MAX;
      for(j=0;j<samples;j++)
      if (sel[j]==1) {
        dis=wdisl1(i,j,sel[j]);
        if (dis<dn) {
          dn=dis;
          jnn=j;
        }
      }

      dd=FLT_MAX;
      for(j=0;j<samples;j++)
      if ((sel[j]==-1)&&(j!=i)) {
        dis=wdisl1(i,j,sel[j]);
        if (dis<dd) {
          dd=dis;
          jdn=j;
        }
      }

      index+=sigmoide(dn/dd);
      expon=exp(rampa-(rampa*(dn/dd)));
      sigmoid=(rampa*expon)/((1+expon)*(1+expon));
      for(k=0;k<dim;k++) {// Update sigmas
        dnk=sigmas[0][k]*fabs(datos[jnn][k]-datos[i][k]);
        ddk=sigmas[1][k]*fabs(datos[jdn][k]-datos[i][k]);
        if ((dnk!=0.0)&&(ddk!=0.0)) {
          // Relevant sigmas
          dif=(datos[jnn][k]-datos[i][k])*(datos[jnn][k]-datos[i][k]);
          fv=sigmas[0][k]*dif/(dnk*dd);
          fv*=sigmoid;
          sg[0][k]+=(mu*fv)/p;
          if (sg[0][k]<0.0) sg[0][k]=0.00001;

          // Non-Relevant sigmas
          dif=(datos[jdn][k]-datos[i][k])*(datos[jdn][k]-datos[i][k]);
          fv=sigmas[1][k]*dif*dn/(ddk*dd*dd);
          fv*=sigmoid;
          sg[1][k]-=(mu*fv)/n;
          if (sg[1][k]<0.0) sg[1][k]=0.00001;
        }
      }
    }
    index/=(p+n);
    copiasig(sigmas,sg);
  }
}
  
double ClassDependentWeightedDistanceQueryCombiner::WeightedL1Distance(float *v1, const vector<double>& v2, int cls) {
  double result=0.0;
  for(int d=0;d<dim;++d) {
    result+=sigmas[cls][d]*fabs(v1[d]-v2[d]);
  }
  return result;
}

void ClassDependentWeightedDistanceQueryCombiner::setParameters(const std::string & parameters) {
    string modestring=getStringAfter(parameters,"MODE=","WDRelevanceScore");
    if(modestring=="WDRelevanceScore") {
      DBG(10) << "mode=WDRelevanceScore" << endl;
      mode_=WDRelevanceScore;
    } else if(modestring=="WDDistSumQuotient") {
      DBG(10) << "mode=WDDistSumQuotient"<<endl;
      mode_=WDDistSumQuotient;
    } else {
      ERR << "Unknown mode. Proceedings with default: WDRelevanceScore" << endl;
      mode_=WDRelevanceScore;
    }
}

float ClassDependentWeightedDistanceQueryCombiner::wdisl1(long int pos1, long int pos2, int sel) {
  float dis;
  int i;
  dis=0;
  int cls;
  if(sel==-1) cls=0;
  else cls=1;
  for (i=0; i<dim; i++) 
    dis+=sigmas[cls][i]*fabs(datos[pos1][i]-datos[pos2][i]);
  return dis;
}

float ClassDependentWeightedDistanceQueryCombiner::sigmoide(float x){
  return 1/(1+exp(rampa-(rampa*x)));
}

string ClassDependentWeightedDistanceQueryCombiner::printSigmas() const {
  ostringstream os;
  for(int j=0;j<2;++j) {
    for(int i=0;i<dim;++i) {
      os << sigmas[j][i] << " ";
    } 
    os <<endl;
  }
  return os.str();
}

RocchioRelevanceFeedbackQueryCombiner::RocchioRelevanceFeedbackQueryCombiner(Retriever& r) : retriever_(r), alpha_(1.0), beta_(-1), gamma_(-1){
}
RocchioRelevanceFeedbackQueryCombiner::~RocchioRelevanceFeedbackQueryCombiner() {
}

void RocchioRelevanceFeedbackQueryCombiner::query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results) {
  ImageContainer Q(*posQ[0]);
  double beta=beta_; double gamma=gamma_;
  
  if(beta<0 and gamma<0) {
    beta=1.0/posQ.size();
    gamma=1.0/negQ.size();  
  }

  
  uint M=Q.numberOfFeatures();
  DBG(10) << "Combining queries into one according to Rocchio" << endl;
  for(uint m=0;m<M;++m) {
    VectorFeature* q=dynamic_cast<VectorFeature*>( Q[m] );
    if(not q) {
      ERR << "Rocchio needs vector features for Q" << endl;
    } else {
      
      for(uint n=0;n<posQ.size();++n) {
        const VectorFeature* v=dynamic_cast<const VectorFeature*>( (*posQ[n])[m] );
        if(not v) {
          ERR << "Rocchio needs vector features for X_p" << endl;
        } else {
          for(uint d=0;d<q->size();++d) {
            (*q)[d]+=beta*(*v)[d];              
          } 
        }
      }

      for(uint n=0;n<negQ.size();++n) {
        const VectorFeature* v=dynamic_cast<const VectorFeature*>( (*negQ[n])[m] );
        if(not v) {
          ERR << "Rocchio needs vector features for X_n" << endl;
        } else {
          for(uint d=0;d<q->size();++d) {
            (*q)[d]-=gamma*(*v)[d];              
          } 
        }
      }      
    } 
  }
  DBG(10) << "Queries combined" << endl;
  
  ScoreSumQueryCombiner adc(retriever_);
  vector<ImageContainer*> positiveQ(1,&Q);
  vector<ImageContainer*> negativeQ(0);
  adc.query(positiveQ, negativeQ, results);
}

void RocchioRelevanceFeedbackQueryCombiner::setParameters(const std::string& parameters) {
  // alpha is currently not being used.
  // default parameters are ok
  alpha_=getDoubleAfter(parameters,"ALPHA=",1.0);
  beta_=getDoubleAfter(parameters,"BETA=",-1);
  gamma_=getDoubleAfter(parameters, "GAMMA=", -1);
}

QueryWeightingQueryCombiner::QueryWeightingQueryCombiner(Retriever& r) : mode_(WeightScores), retriever_(r), FLT_MAX(100000000.0), FLT_MIN(-100000000.0), rampa(10.0), posWeight_(1.0), negWeight_(0.8333){
}

QueryWeightingQueryCombiner::~QueryWeightingQueryCombiner() {
}


double QueryWeightingQueryCombiner::distBetweenImages(const ImageContainer* q1, const ImageContainer* q2) {  
  return log(retriever_.scorer()->getScore(retriever_.imageComparator().compare(q1, q2)));
}

void QueryWeightingQueryCombiner::calcWeights(const std::vector<ImageContainer*>& posQ, const std::vector<ImageContainer*>& negQ, vector<double>& posWeights, vector<double>& negWeights) {
  int it,maxit=100;
  int jdn=0,jnn=0;
  float index,dn,dd,fv,dis;
  float mu=.1,expon,sigmoid;
  int p=posQ.size(),n=negQ.size();
  
  if ((n==0)||(p<2)) return;

  vector<double> posSg(posWeights), negSg(negWeights);
  DBG(50) << "PosW";    for (uint i=0; i<posSg.size(); ++i) {      BLINK(50) <<" "<< posSg[i];    }    BLINK(50) << endl;
  DBG(50) << "negW";    for (uint i=0; i<negSg.size(); ++i) {      BLINK(50) <<" "<< negSg[i];    }    BLINK(50) << endl;

  it=0;
  while (it<maxit) {
    DBG(50) << VAR(it) << endl;
    vector<double> posSg(posWeights), negSg(negWeights);  
    DBG(50) << "PosW";    for (uint i=0; i<posSg.size(); ++i) {      BLINK(50) <<" "<< posSg[i];    }    BLINK(50) << endl;
    DBG(50) << "negW";    for (uint i=0; i<negSg.size(); ++i) {      BLINK(50) <<" "<< negSg[i];    }    BLINK(50) << endl;

    it++;
    index=0.0;

    for(int i=0;i<p;++i) {
      dn=FLT_MAX;
      for (int j=0; j<p; j++) {
        if (j!=i) {
          dis=distBetweenImages(posQ[i], posQ[j]);
          if (dis<dn) {
            dn=dis;
            jnn=j;
          }
        }
      }

      dd=FLT_MAX;
      for(int j=0;j<n;j++) {
        dis=distBetweenImages(posQ[i],negQ[j]);
        if (dis<dd) {
          dd=dis;
          jdn=j;
        }
      }
      index+=sigmoide(dn/dd)/p;
      expon=exp(rampa-(rampa*(dn/dd)));
      sigmoid=(rampa*expon)/((1+expon)*(1+expon));

      fv=(sigmoid*dn/dd)/posWeights[jnn];
      posSg[jnn]-=fv*mu/p;
      if (posSg[jnn]<0.0) posSg[jnn]=0.0000001;

      fv=(sigmoid*dn/dd)/negWeights[jdn];
      negSg[jdn]+=fv*mu/p;
    }


    for (int i=0; i<n; ++i) {
      dn=FLT_MAX;
      for (int j=0; j<p; j++) {
        dis=distBetweenImages(negQ[i], posQ[j]);
        if (dis<dn) {
          dn=dis;
          jnn=j;
        }
      }
    

    dd=FLT_MAX;
    for(int j=0;j<n;j++) {        
      if (j!=i) {
        dis=distBetweenImages(negQ[i],negQ[j]);
        if (dis<dd) {
          dd=dis;
          jdn=j;
        }
      }
    }

    index-=sigmoide(dn/dd)/n;
    expon=exp(rampa-(rampa*(dn/dd)));
    sigmoid=(rampa*expon)/((1+expon)*(1+expon));

    fv=(sigmoid*dn/dd)/posWeights[jnn];
    posSg[jnn]+=fv*mu/n;

    fv=(sigmoid*dn/dd)/negWeights[jdn];
    negSg[jdn]-=fv*mu/n;
    if (negSg[jdn]<0.0) negSg[jdn]=0.0000001;
    }
    posWeights=posSg; negWeights=negSg;
  }
}

void QueryWeightingQueryCombiner::query(const std::vector<ImageContainer*> posQueries, const std::vector<ImageContainer*> negQueries, std::vector<ResultPair>& results) {
  vector<double> posWeights(posQueries.size(),1.0),negWeights(negQueries.size(),1.0);
  double posWeight=posWeight_/posQueries.size();
  double negWeight=negWeight_/negQueries.size();

  retriever_.imageComparator().tuneDistances(posQueries,negQueries);
  calcWeights(posQueries, negQueries, posWeights, negWeights);
  
  uint N=retriever_.database().size();

  vector<double> activeScores(N, 0.0);
  DBG(10) << "Querying ";
  for (long q=0; q<long(posQueries.size()); ++q) {
    BLINK(10) << "+";
    DBG(15) << "Positive query: " << posQueries[q]->basename() << endl;
    if (retriever_.imageComparator().size() != posQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than positive query " << q << " (" << posQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(posQueries[q]);
    retriever_.getScores(posQueries[q], activeScores);
    retriever_.imageComparator().stop();
    if(mode_==WeightScores) {
      for (uint i=0; i<N; ++i) {
        results[i].first+=posWeights[q]*posWeight*activeScores[i];
      }
    } else if(mode_==WeightDistances) {
      for (uint i=0; i<N; ++i) {
        results[i].first+=posWeights[q]*exp(posWeight*log(activeScores[i]));
      }    
    }
  }

  for (long q=0; q<long(negQueries.size()); ++q) {
    BLINK(10) << "-";
    DBG(15) << "Negative query: " << negQueries[q]->basename() << endl;

    if (retriever_.imageComparator().size() != negQueries[q]->numberOfFeatures()) {
      ERR << "ImageComparator has different number of distances (" << retriever_.imageComparator().size() << ") than negative query " << q << " (" << negQueries[q]->numberOfFeatures() << ")." << endl;
    }

    retriever_.imageComparator().start(negQueries[q]);
    retriever_.getScores(negQueries[q], activeScores);
    retriever_.imageComparator().stop();
    if(mode_==WeightScores) {
      for (uint i=0; i<N; ++i) {
        results[i].first+=negWeights[q]*negWeight*(1.0-activeScores[i]);
      }
    } else if(mode_==WeightDistances) {
      for (uint i=0; i<N; ++i) {
        results[i].first+=negWeight*(1.0-log(exp(negWeights[q]*activeScores[i])));
      }
    }
  }
  BLINK(10) << endl;
}

float QueryWeightingQueryCombiner::sigmoide(float x){
  return 1/(1+exp(rampa-(rampa*x)));
}


void QueryWeightingQueryCombiner::setParameters(const std::string& parameters) {
  string modestring=getStringAfter(parameters,"MODE=","WeightScores");
  if(modestring=="WeightScores") {DBG(10) << "mode=WeightScores" << endl; mode_=WeightScores;}
  else if(modestring=="WeightDistances") {DBG(10) << "mode=WeightDistances" << endl; mode_=WeightDistances;}
  else {ERR << "Unknown Mode. Using default: WeightScores" << endl; mode_=WeightScores;}
}

  
