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
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stack>
#include <algorithm>
#include "retriever.hpp"
#include "dist_metafeature.hpp"
#include "dist_textfeature.hpp"
#include "textfeature.hpp"
#include "getscoring.hpp"
#include "net.hpp"


using namespace std;

#include "ScopeTimer.h"

Retriever::Retriever() :
  database_(), imageComparator_(), queryCombiner_(new ScoreSumQueryCombiner(*this)), reRanker_(new ReRanker(*this)), results_(0), extensions_(0), interactor_(), filterApply_(false), partialLoadingApply_(false), filter_() {
  scorer_=new LinearScoring();
  //  queryCombiner_=new AddingQueryCombiner(*this);
}

Retriever::~Retriever() {
  delete scorer_;
  delete queryCombiner_;
  database_.clear();
}

void Retriever::initialize() {
  imageComparator_.initialize(database_);
}

void Retriever::setScoring(const string &scoringname) {
  delete scorer_;
  scorer_=getScoring(scoringname, database_.numberOfSuffices());
}

void Retriever::setQueryCombiner(const string &queryCombiningName) {
  delete queryCombiner_;
  if (queryCombiningName.find("scoresum")==0) {
    DBG(10) << "AddingQueryCombiner" << endl;
    queryCombiner_=new ScoreSumQueryCombiner(*this);
  } else if (queryCombiningName.find("nnquot")==0) {
    queryCombiner_=new NNQuotientQueryCombiner(*this);
  } else if (queryCombiningName.find("sumquot")==0) {
    queryCombiner_=new SumQuotientQueryCombiner(*this);
  } else if (queryCombiningName.find("relevanceScore")==0) {
    queryCombiner_=new RelevanceScoreQueryCombiner(*this);
  } else if (queryCombiningName.find("distsumquot")==0) {
    queryCombiner_=new DistSumQuotientQueryCombiner(*this);
  } else if (queryCombiningName.find("svm")==0) {
    queryCombiner_=new SVMQueryCombiner(*this);
  } else if (queryCombiningName.find("wd")==0) {
    queryCombiner_=new WeightedDistanceQueryCombiner(*this);
  } else if (queryCombiningName.find("cdwd")==0) {
    queryCombiner_=new ClassDependentWeightedDistanceQueryCombiner(*this);
  } else if (queryCombiningName.find("rocchio")==0) {
    queryCombiner_=new RocchioRelevanceFeedbackQueryCombiner(*this);
  } else if (queryCombiningName.find("queryweighting")==0) {
    queryCombiner_=new QueryWeightingQueryCombiner(*this);
  } else {
    ERR << "Unknown query combiner '"<<queryCombiningName<<"'. Using 'adding' with default parameters." <<endl;
    queryCombiner_=new ScoreSumQueryCombiner(*this);
  }
  queryCombiner_->setParameters(queryCombiningName);
}

void Retriever::setReranking(const string &rerankingName){

  delete reRanker_;
  if (rerankingName.find("none")==0){
    reRanker_=new ReRanker(*this);
  } else if (rerankingName.find("cluster")==0){
      DBG(10)<<"Adding clustering reranker " << endl;
      reRanker_=new ClusterReranker(*this);
  } else if (rerankingName.find("greedy")==0){
    DBG(10)<<"Adding greedy reranker "<<endl;
    reRanker_=new GreedyReranker(*this);
  } else if (rerankingName.find("dp")==0) {
    DBG(10) << "DPOptimisingReranker" << endl;
    reRanker_=new DPOptimisingReranker(*this);
  } else if (rerankingName.find("monoDP")==0) {
    DBG(10) << "DP2NeyReranker" << endl;
    reRanker_=new DPOpt2(*this);
  }  else {
    ERR << "Unknown reranking '"<<rerankingName<<"'. Using 'none' with default parameters."<<endl;
    reRanker_=new ReRanker(*this);
  }
  reRanker_->setParameters(rerankingName);
}
  

void Retriever::resolveNames(const vector< string > &queryNames, vector<ImageContainer*> &queries, stack<ImageContainer*> &newCreated) {
  for (uint i=0; i<queryNames.size(); ++i) {
    ImageContainer* tmp=NULL;
    tmp=database_.getByName(queryNames[i]);
    if (!tmp) {
      DBG(10) << "Image not in database, trying to load: " << queryNames[i];
      tmp=new ImageContainer(queryNames[i],database_.numberOfSuffices());
      if (!database_.loadQuery(queryNames[i], tmp)) {
        delete tmp;
        tmp=NULL;
      }
      if (tmp) {
        newCreated.push(tmp);
        BLINK(10) << " ... done"<< endl;
      } else
        BLINK(10) << "... failed" << endl;
    }
    if (tmp)
      queries.push_back(tmp);
  }
}

void Retriever::retrieve(const vector< string >& posQueryNames, const vector< string >& negQueryNames, vector<ResultPair>& results) {

  // get image containers for these images
  vector<ImageContainer*> posQueries;
  vector<ImageContainer*> negQueries;

  stack<ImageContainer*> newCreated; // memorize those images that must be deleted after processing the query

  resolveNames(posQueryNames, posQueries, newCreated);
  resolveNames(negQueryNames, negQueries, newCreated);

  retrieve(posQueries, negQueries, results);
  
  vector<ResultPair> tmp;
  reRanker_->rerank(posQueries, negQueries, results,tmp);
  results=tmp;
  
  while (!newCreated.empty()) {
    delete newCreated.top();
    newCreated.pop();
  }
}

void Retriever::getScores(const ImageContainer* q, vector<double>&scores) {

  ScopeTimer st1((char*)"Retriever::getScores");

  uint N=database_.size();
  uint M=database_.numberOfSuffices();

  vector< vector<double> > distMatrix(N, vector<double>(M));
  vector<double> imgDists;

  // from here, we want parallelization using OpenMP
#pragma omp parallel
  {
    //get distance to each of the database images
    { // begin "get distance" scope
      ScopeTimer st2((char*)"Retriever::getScores -> get distance");

      //this next llop can be done in parallel
      //but each of the threads needs its own imgDists vector
      //all the other variables (q,database_[i]) are readonly
      //when doing things parallel: take care with static, shared mem
#pragma omp for schedule(static) private(imgDists)
      for (long i=0; i<long(N); ++i) {
        vector<double>&d=distMatrix[i];
        imgDists=imageComparator_.compare(q, database_[i]);
        for (long j=0; j<long(M); ++j) {
          d[j]=imgDists[j];
        }
      }

      //here we wait until all threads have finished.
      //that is, distMatrix is completely
      //this is called barrier; to enforce this manually: #pragma omp barrier
    } // end "get distance" scope


    //normalize the distances
    { // begin "normalize" scope
      ScopeTimer st3((char*)"Retriever::getScores -> normalize");

      //this next loop is parallelized, too. note that loop-local
      //variables are not problematic, as they are instantiated for each
      //thread
#pragma omp for schedule(static)
      for (long j=0; j<long(M); ++j) {
        double sum=0.0;
        for (long i=0; i<long(N); ++i) {
          sum+=distMatrix[i][j];
        }
        sum/=double(N);
        if (sum!=0.0) {
          double tmp=1/sum;
          for (long i=0; i<long(N); ++i) {
            distMatrix[i][j]*=tmp;
          }
        }
      }
    } // end "normalize" scope
  } // end omp parallel

  // here: distance interactions: this is still quite buggy
  { // begin "interactor" scope
    ScopeTimer st4((char*)"Retriever::getScores -> interactor.apply()");
    interactor_.apply(distMatrix);
  } // end "interactor" scope

  //now get the scores
  { // begin "get the scores" scope
    ScopeTimer st5((char*)"Retriever::getScores -> get the scores");
    for (uint i=0; i<N; ++i) {
      DBG(105) << "DISTS:";
      for (uint j=0; j<M; ++j) {
        BLINK(105) << distMatrix[i][j] << " ";
      }
      scores[i]=scorer_->getScore(distMatrix[i]);
      BLINK(105) << "->" << scores[i] << endl;
    }
  } // end "get the scores" scope
}

void Retriever::saveDistances(string imagename, string filename) {
  /*----------------------------------------------------------------------
   * get distance matrix 
   * --------------------------------------------------------------------*/
  uint N=database_.size();
  uint M=database_.numberOfSuffices();
  bool newlyLoaded=false;
  vector< vector<double> > distMatrix(N, vector<double>(M));
  vector<double> imgDists;

  ImageContainer *q=database_.getByName(imagename);
  if (!q) {
    DBG(10) << "Image not in database, trying to load: " << imagename;
    q=new ImageContainer(imagename,database_.numberOfSuffices());
    if (!database_.loadQuery(imagename, q)) {
      delete q;
      q=NULL;
      BLINK(10) << "... failed" << endl;
    } else {
      BLINK(10) << " ... done"<< endl;
      newlyLoaded=true;
    }
  }
  imageComparator_.start(q);
  //get distance to each of the database images
  for (uint i=0; i<N; ++i) {
    vector<double>&d=distMatrix[i];
    imgDists=imageComparator_.compare(q, database_[i]);
    for (uint j=0; j<M; ++j) {
      d[j]=imgDists[j];
    }
  }
  imageComparator_.stop();

  //normalize
  for (uint j=0; j<M; ++j) {
    double sum=0.0;
    for (uint i=0; i<N; ++i) {
      sum+=distMatrix[i][j];
    }
    sum/=double(N);
    if (sum!=0.0) {
      double tmp=1/sum;
      for (uint i=0; i<N; ++i) {
        distMatrix[i][j]*=tmp;
      }
    }
  }

  /*----------------------------------------------------------------------
   * save distance matrix
   * --------------------------------------------------------------------*/

  ifstream is;
  is.open(filename.c_str());
  if (is.good()) {
    ERR << "Distance file to be written '" << filename << "' exists already. Skipping." << endl;
  } else {

    ofstream os;
    os.open(filename.c_str());
    if (!os.good()) {
      ERR << "Cannot write distance file '" <<filename << "'." << endl;
      exit(10);
    } else {
      os << "# distmatrix for "<< q->basename() << endl;
      os << "# distances";
      for (uint j=0; j<M; ++j) {
        os << " " << imageComparator_.distance(j)->name();
      }
      os << endl;
      os << "nofdistances "<< N << " " << M << endl;

      for (uint i=0; i<N; ++i) {
        os << i;
        for (uint j=0; j<M; ++j) {
          os << " "<<distMatrix[i][j];
        }
        os << endl;
      }
      os.close();
    }

    if (newlyLoaded) {
      delete q;
      q=NULL;
    }
  }
}

// fill the distance matrix taking into account that some distance have NOT been calculated.
// put 1.2*maxdist into these positions
void Retriever::getDistances(const ImageContainer* q, const vector<uint>& stillToConsider, const vector<uint>& depreciated, vector< vector<double> >& distMatrix, const uint distanceID) {
  double imgDist;
  double maxDist = 0.0;

  for (long i=0; i<long(stillToConsider.size()); ++i) {
    imgDist = imageComparator_.compare(q, database_[stillToConsider[i]], distanceID);
    distMatrix[stillToConsider[i]][distanceID] = imgDist;
    maxDist = max(maxDist, imgDist);
  }

  // beware: here an overflow might occur !
  maxDist *= 1.2;

  // now we set the distance matrix values for the unwanted images
  for (long j=0; j<long(depreciated.size()); ++j) {
    distMatrix[depreciated[j]][distanceID] = maxDist;
  }
}

void Retriever::getBest(vector<uint> &stillToConsider, vector<uint> &depreciated, const vector<double> &scores, uint &amount) {

  uint imageSize = stillToConsider.size();
  vector<ResultPair> toSort;
  uint img;
  for (uint i=0; i<imageSize; ++i) {
    img = stillToConsider[i];
    toSort.push_back(ResultPair(scores[img], img));
  }

  // sorting indizies according to scores
  sort(toSort.rbegin(), toSort.rend());

  if (amount > imageSize) {
    amount = imageSize;
  }

  // tmp construct
  vector<uint> imageTmp;
  for (uint j=0; j<amount; ++j) {
    img = toSort[j].second;
    imageTmp.push_back(img);
  }

  for (uint i=amount; i<imageSize; ++i) {
    img = toSort[i].second;
    depreciated.push_back(img);
  }

  // swap contents of tmp construct and real construct
  swap(stillToConsider, imageTmp);
}

void Retriever::getScores(vector<vector<double> > &distMatrix, vector<double> &scores) {
  uint N=database_.size();
  uint M=database_.numberOfSuffices();

  scores=vector<double>(N, 0.0);

  for (long j=0; j<long(M); ++j) {
    double sum=0.0;
    for (long i=0; i<long(N); ++i) {
      sum+=distMatrix[i][j];
    }
    sum/=double(N);
    if (sum!=0.0) {
      double tmp=1/sum;
      for (long i=0; i<long(N); ++i) {
        distMatrix[i][j]*=tmp;
      }
    }
  }

  interactor_.apply(distMatrix);

  for (uint i=0; i<N; ++i) {
    DBG(105) << "DISTS:";
    for (uint j=0; j<M; ++j) {
      BLINK(105) << distMatrix[i][j] << " ";
    }
    scores[i]=scorer_->getScore(distMatrix[i]);
    BLINK(105) << "->" << scores[i] << endl;
  }
}

void Retriever::retrieve(const vector<ImageContainer*>& posQueries, const vector<ImageContainer*>& negQueries, vector<ResultPair>& results) {

  uint N=database_.size();
  vector<double> activeScores(N, 0.0);
  //init result field
  results.clear();
  for (uint i=0; i<database_.size(); ++i) {
    results.push_back(ResultPair(0.0, i));
  }
  if (!filterApply_) {
    //positive queries

    queryCombiner_->query(posQueries, negQueries, results);

    //check whether query expansion has to be done
    if (extensions_!=0) {
      //if query expansion has to be done:
      // get ranks
      sort(results.rbegin(), results.rend());

      vector<ImageContainer*> expansion;

      //copy extensions_ into positive queries
      for (uint i=0; i<extensions_; ++i) {
        expansion.push_back(database_[results[i].second]);
      }

      //init result field
      results.clear();
      for (uint i=0; i<N; ++i) {
        results.push_back(ResultPair(0.0, i));
      }

      // and requery using these positive queries
      for (uint q=0; q<expansion.size(); ++q) {
        imageComparator_.start(expansion[q]);
        getScores(expansion[q], activeScores);
        imageComparator_.stop();
        for (uint i=0; i<N; ++i) {
          results[i].first+=activeScores[i];
        }
      }
    }
  } else { // filterApply=true

    // the stillToConsider vector contains the best images as results of the filtered retrieval
    // and this vector is reused through out the filtered retrieval
    // in analogous fashion the neqimages vector contains the remaining database images
    vector<uint> stillToConsider;
    vector<uint> depreciated;

    // initial value for the distances. this value is used to ensure that the distance matrix
    // contains in every step definied values and not just memory garbage
    const double initDummyDist = 100.0;
    uint lbffidx=0;
    uint M = database_.numberOfSuffices();

    //positive queries
    for (long q=0; q<long(posQueries.size()); ++q) {

      vector< vector <double> > distMatrix(N, vector<double>(M, initDummyDist));
      for (uint i=0; i<N; ++i) {
        stillToConsider.push_back(i);
      }
      vector<double> activeScores(N, 0.0);

      DBG(10) << "Positive query: " << posQueries[q]->basename() << endl;
      if (imageComparator_.size() != posQueries[q]->numberOfFeatures()) {
        ERR << "ImageComparator has different number of distances (" << imageComparator_.size() << ") than positive query " << q << " (" << posQueries[q]->numberOfFeatures() << ")." << endl;
      }
      for (uint i=0; i<filter_.size(); ++i) {
        DBG(20) << "Filtering " << i+1 << "/" << filter_.size()<< ":  " <<filter_[i].first << ":" << filter_[i].second << endl;
        // first check whether or not feature information has to be loaded into
        // the database
        if (partialLoadingApply_) {
          // check if the current feature is not loaded
          lbffidx=filter_[i].first;
          if (database_.binFilesNotToLoad(lbffidx)) {
            DBG(105) << "DEBUG Load from LBFF with idx " << lbffidx << endl;
            database_.loadFromLBFF(lbffidx, stillToConsider);
          }
        }
        imageComparator_.start(posQueries[q], filter_[i].first);
        getDistances(posQueries[q], stillToConsider, depreciated, distMatrix, filter_[i].first);
        imageComparator_.stop(filter_[i].first);
        getScores(distMatrix, activeScores);
        // remove the loaded feature information if partial loading
        // doing this as early as possible
        if (partialLoadingApply_ && database_.binFilesNotToLoad(lbffidx)) {
          database_.removeFeatureInformation(lbffidx, stillToConsider);
        }
        getBest(stillToConsider, depreciated, activeScores, filter_[i].second);
      }

      for (uint i=0; i<N; ++i) {
        results[i].first+=activeScores[i];
      }

      // cleaning memory for next iteration
      stillToConsider.clear();
      depreciated.clear();
      distMatrix.clear();
    } // end for loop for positive queries

    //negative queries
    for (long q=0; q<long(negQueries.size()); ++q) {

      vector< vector <double> > distMatrix(N, vector<double>(M, initDummyDist));
      for (uint i=0; i<N; ++i) {
        stillToConsider.push_back(i);
      }
      vector<double> activeScores(N, 0.0);

      DBG(10) << "Negative query: " << negQueries[q]->basename() << endl;
      if (imageComparator_.size() != negQueries[q]->numberOfFeatures()) {
        ERR << "ImageComparator has different number of distances (" << imageComparator_.size() << ") than negative query " << q << " (" << negQueries[q]->numberOfFeatures() << ")." << endl;
      }
      for (uint i=0; i<filter_.size(); ++i) {
        // first check whether or not feature information has to be loaded into
        // the database
        if (partialLoadingApply_) {
          // check if the current feature is not loaded
          lbffidx=filter_[i].first;
          if (database_.binFilesNotToLoad(lbffidx)) {
            database_.loadFromLBFF(lbffidx, stillToConsider);
          }
        }
        imageComparator_.start(negQueries[q], filter_[i].first);
        getDistances(negQueries[q], stillToConsider, depreciated, distMatrix, filter_[i].first);
        imageComparator_.stop(filter_[i].first);
        getScores(distMatrix, activeScores);
        if (partialLoadingApply_ && database_.binFilesNotToLoad(lbffidx)) {
          database_.removeFeatureInformation(lbffidx, stillToConsider);
        }
        getBest(stillToConsider, depreciated, activeScores, filter_[i].second);
      }

      for (uint i=0; i<N; ++i) {
        results[i].first+=(1-activeScores[i]);
      }

      // cleaning memory for next iteration
      stillToConsider.clear();
      depreciated.clear();
      distMatrix.clear();
    } // end negative query

    //check whether query expansion has to be done
    if (extensions_!=0) {
      //if query expansion has to be done:
      // get ranks
      sort(results.rbegin(), results.rend());

      vector<ImageContainer*> expansion;

      //copy extensions_ into positive queries
      for (uint i=0; i<extensions_; ++i) {
        expansion.push_back(database_[results[i].second]);
      }

      //init result field
      results.clear();
      for (uint i=0; i<N; ++i) {
        results.push_back(ResultPair(0.0, i));
      }

      // and requery using these positive queries
      for (uint q=0; q<expansion.size(); ++q) {

        vector< vector <double> > distMatrix(N, vector<double>(M, initDummyDist));
        for (uint i=0; i<N; ++i) {
          stillToConsider.push_back(i);
        }

        for (uint i=0; i<filter_.size(); ++i) {
          // first check whether or not feature information has to be loaded into
          // the database
          if (partialLoadingApply_) {
            // check if the current feature is not loaded
            lbffidx=filter_[i].first;
            if (database_.binFilesNotToLoad(lbffidx)) {
              database_.loadFromLBFF(lbffidx, stillToConsider);
            }
          }
          imageComparator_.start(expansion[q], filter_[i].first);
          getDistances(expansion[q], stillToConsider, depreciated, distMatrix, filter_[i].first);
          imageComparator_.stop(filter_[i].first);
          getScores(distMatrix, activeScores);
          if (partialLoadingApply_ && database_.binFilesNotToLoad(lbffidx)) {
            database_.removeFeatureInformation(lbffidx, stillToConsider);
          }
          getBest(stillToConsider, depreciated, activeScores, filter_[i].second);
        }

        for (uint i=0; i<N; ++i) {
          results[i].first+=activeScores[i];
        }

        // cleaning memory for next iteration
        stillToConsider.clear();
        depreciated.clear();
        distMatrix.clear();
      }
    } // end extensions
  } // end else
  DBG(15) << "end retrieve" << endl;
}

vector<ResultPair> Retriever::metaretrieve(const string& query) {
  uint N=database_.size();
  uint M=database_.numberOfSuffices();
  vector<ResultPair> result;

  //init result field
  /*  for(uint i=0;i<database_.size();++i) {
   result.push_back(ResultPair(0.0,i));
   }*/

  MetaFeatureDistance metadist;

  // Find index of metafeature
  string dist_name;
  int metafeatureidx = -1;
  for (unsigned int i=0; i<M; ++i) {
    dist_name = imageComparator_.distance(i)->name();
    if (dist_name=="metafeature") {
      metafeatureidx = i;
      break;
    }
  }

  if (metafeatureidx==-1) {
    return vector<ResultPair>(0);
  }

  //Now parse the query string and build the corresponding metafeature
  ImageContainer imgcon("temporary query object", M);

  vector<string> tokens;
  const string delimiters = ":,";

  // Tokenizer ripped from the web (what a shame that c++ doesn't have one)
  // Skip delimiters at beginning.
  string::size_type lastPos = query.find_first_not_of(delimiters, 0);

  // Find first "non-delimiter".
  string::size_type pos = query.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(query.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = query.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = query.find_first_of(delimiters, lastPos);
  }

  // Quick sanity check: even number of tokens
  if (tokens.size() % 2 > 0) {
    DBG(10) << "Invalid metaquery!" << endl;
    return result;
  }

  // Another sanity check: only valid keys in the query
  ::std::vector< ::std::string > valid_keys;
  valid_keys = getMetaFeatureInfo().first;
  for (unsigned i=0; i<tokens.size(); i+=2) {
    string key = tokens[i];
    if (find(valid_keys.begin(), valid_keys.end(), key) == valid_keys.end()) {
      // Alert!
      DBG(10) << "Invalid metaquery!" << endl;
      return result;
    }
  }

  // Now make the metafeature map from the token list
  std::map<std::string,std::string> val;

  string dbgstr = "Metaquery - ";
  string key, value;
  vector<string>::iterator t_it;
  for (t_it = tokens.begin(); t_it!=tokens.end(); ++t_it) {
    key = string(*t_it);
    ++t_it;
    value = string(*t_it);
    val[key] = value;
    dbgstr += key + ":" + value + " ";
  }
  DBG(10) << dbgstr << endl;

  // This is a virtual image that only has a meta feature
  imgcon[metafeatureidx] = new MetaFeature(val);

  vector<double> distsToImages(N);
  vector<double> imgDists;

  // get distance to each of the database images
  for (uint i=0; i<N; ++i) {
    BaseFeature *dbmf= database_[i]->operator[](metafeatureidx);
    distsToImages[i]=metadist.distance(imgcon[metafeatureidx], dbmf);
  }

  // normalize
  double sum=0.0;
  for (uint i=0; i<N; ++i) {
    sum+=distsToImages[i];
  }
  sum/=double(N);
  if (sum!=0.0) {
    double tmp=1/sum;
    for (uint i=0; i<N; ++i) {
      distsToImages[i]*=tmp;
    }
  }

  // Compute the scores (this omits the use of a scoring class)
  unsigned max_dist = val.size();
  for (uint i=0; i<N; ++i) {
    if (distsToImages[i] < max_dist) {
      result.push_back(ResultPair(exp(-distsToImages[i]), i));
    }
  }

  //check whether query expansion has to be done
  if (extensions_!=0) {
    vector<double> activeScores(N, 0.0);
    //if query expansion has to be done:

    // get ranks
    sort(result.rbegin(), result.rend());

    //reset positive queries (also clear negative queries)
    vector<ImageContainer*> expansion;

    //copy extensions_ into positive queries
    for (uint i=0; i<extensions_; ++i) {
      expansion.push_back(database_[result[i].second]);
    }

    //init result field
    result.clear();
    for (uint i=0; i<database_.size(); ++i) {
      result.push_back(ResultPair(0.0, i));
    }

    // and requery using these positive queries
    for (uint q=0; q<expansion.size(); ++q) {
      imageComparator_.start(expansion[q]);
      getScores(expansion[q], activeScores);
      imageComparator_.stop();
      for (uint i=0; i<N; ++i) {
        result[i].first+=activeScores[i];
      }
    }
  }

  return result;
}

vector<ResultPair> Retriever::textretrieve(const string& query) {
  uint N=database_.size();
  uint M=database_.numberOfSuffices();
  vector<ResultPair> result;

  // Where is which textfeature and which langage do they have ?
  string dist_name;
  TextFeatureDistance *tfd;
  map<string,TextFeatureDistance*> textdists;
  map<string,unsigned> textdistindices;
  for (unsigned int i=0; i<M; ++i) {
    dist_name = imageComparator_.distance(i)->name();
    if (dist_name=="textfeature") {
      tfd = (TextFeatureDistance*)imageComparator_.distance(i);
      textdists[tfd->language()] = tfd;
      textdistindices[tfd->language()] = i;
      DBG(10) << "found " << tfd->language() << " at index " << i << endl;
    }
  }

  // Decompose query string in its languages
  DBG(15) << "decomposing query" << endl;
  map<string,string> cl_queries;
  map<string,TextFeatureDistance*>::iterator tdi;
  for (tdi=textdists.begin(); tdi!=textdists.end(); ++tdi) {
    string lang = tdi->first;
    string langpart;
    unsigned langpos = query.find(lang);
    unsigned langbegin=0, langend=0;
    if (langpos != string::npos) {
      // Sanity check
      if (query.substr(langpos + lang.size(), 2) != ":\"") {
        DBG(10) << "Malformed query for " << lang << ": " << query << endl;
        return vector<ResultPair>(0);
      }

      langbegin = langpos + lang.size() + 2;
      langend = query.find('"', langbegin);
      if (langend == string::npos) {
        DBG(10) << "Missing closing quotation mark for " << lang << endl;
        return vector<ResultPair>(0);
      }
      langpart = query.substr(langbegin, langend - langbegin);
      cl_queries[lang] = langpart;
      DBG(10) << "Query in " << lang << ": " << langpart << endl;
    }
  }

  // Sanity check: Are there no cl-queries although there are
  // specified languages for the textfeature(s)?
  if ((cl_queries.size()==0)&&(textdists.begin()->first!="None")) {
    DBG(10) << "Error: User did not specify language!" << endl;
    return vector<ResultPair>(0);
  }

  // If nothing is found, but the string is non-empty, it's a
  // normal query, not a cross-language one.
  if (cl_queries.size()==0) {
    DBG(10) << "normal query" << endl;
    TextFeatureDistance *textdist;
    textdist = textdists["None"];

    // Do the textretriever query
    textdist->query_wmir(query);

    // Find index of textfeature
    // (we assume there is only one, but if there are multiple ones,
    // the last one is taken (for no particular reason))
    int textfeatureidx = -1;
    for (unsigned int i=0; i<M; ++i) {
      dist_name = imageComparator_.distance(i)->name();
      if (dist_name=="textfeature") {
        textfeatureidx = i;
        break;
      }
    }

    // Make an empty image
    ImageContainer imgcon("temporary query object", M);
    imgcon[textfeatureidx] = new TextFeature();

    vector<double> distsToImages(N);

    // get distance to each of the database images
    BaseFeature *dbmf;

    // Find the maximum occuring distance in textretrieve. Assume that
    // these images have not been returned by the text retrieval
    // engine. This value is necessary for discarding these images later.
    double maxdist_emp = 0.0;

    for (uint i=0; i<N; ++i) {
      dbmf = database_[i]->operator[](textfeatureidx);
      distsToImages[i]=textdist->distance(imgcon[textfeatureidx], dbmf);
      if (maxdist_emp<distsToImages[i]) {
        maxdist_emp=distsToImages[i];
      }
    }

    // normalize
    double sum=0.0;
    double normalization_factor = 1.0;
    for (uint i=0; i<N; ++i) {
      sum+=distsToImages[i];
    }
    sum/=double(N);
    if (sum!=0.0) {
      normalization_factor = 1/sum;
      for (uint i=0; i<N; ++i) {
        distsToImages[i]*=normalization_factor;
      }
    }

    // Compute the scores (this omits the use of a scoring class)

    // maxval is the normalized dist value for an image that was not
    // in the textretriever's result.
    double maxval = maxdist_emp * normalization_factor;

    for (uint i=0; i<N; ++i) {
      if (distsToImages[i] != maxval) {
        ResultPair r;
        r.first = exp(-distsToImages[i]);
        r.second = i;
        result.push_back(r);
      }
    }

  } else {
    cout << "x-language query" << endl;

    // Query all distances with the correct queries
    map<string,string>::iterator clqi;
    for (clqi=cl_queries.begin(); clqi!=cl_queries.end(); ++clqi) {
      cout << "querying " << clqi->first << endl;
      textdists[clqi->first]->query_wmir(clqi->second);
    }

    // Make new maps that only contain the textdists and indices of languages that were
    // used for querying
    map<string,TextFeatureDistance*> querytextdists;
    map<string,unsigned> querytextdistindices;
    map<string,TextFeatureDistance*>::iterator tdi;
    map<string,unsigned>::iterator tdii;
    cout << "Only calculating dists for" << endl;
    for (tdi=textdists.begin(), tdii=textdistindices.begin(); tdi != textdists.end(); ++tdi, ++tdii) {
      clqi = cl_queries.find(tdi->first); // Find language in the query
      if (clqi != cl_queries.end()) { // If this language was part of the query, add it
        querytextdists[tdi->first] = tdi->second;
        querytextdistindices[tdii->first] = tdii->second;
        cout << tdii->first << endl;
      }
    }

    // Make an empty image with the queried textfeatures
    DBG(10) << "Making empty image" << endl;
    ImageContainer imgcon("temporary query object", M);
    for (tdii=querytextdistindices.begin(); tdii!=querytextdistindices.end(); ++tdii) {
      imgcon[tdii->second] = new TextFeature();
    }

    // get distance to each of the database images
    DBG(10) << "Getting distances" << endl;
    BaseFeature *dbmf;
    map<string, vector<double> > distsToImages;
    for (tdii=querytextdistindices.begin(); tdii!=querytextdistindices.end(); ++tdii) {
      string tdlang = tdii->first;
      unsigned tdidx = tdii->second;

      distsToImages[tdlang].resize(N, 0.0);
      for (uint i=0; i<N; ++i) {
        dbmf = database_[i]->operator[](tdidx);
        distsToImages[tdlang][i]=querytextdists[tdlang]->distance(imgcon[tdidx], dbmf);
      }
    }

    // normalize
    DBG(10) << "normalizing" << endl;
    map<string,double> normalization_factors;
    for (tdii=querytextdistindices.begin(); tdii!=querytextdistindices.end(); ++tdii) {
      double sum=0.0;
      normalization_factors[tdii->first] = 1.0;

      for (uint i=0; i<N; ++i) {
        sum+=distsToImages[tdii->first][i];
      }
      sum/=double(N);
      if (sum!=0.0) {
        normalization_factors[tdii->first] = 1/sum;
        for (uint i=0; i<N; ++i) {
          distsToImages[tdii->first][i]*=normalization_factors[tdii->first];
        }
      }
    }

    // Compute the scores (this omits the use of a scoring class)

    // maxval is the normalized dist value for an image that was not
    // in the textretriever's result.
    map<string, double> maxvals;
    for (tdii=querytextdistindices.begin(); tdii!=querytextdistindices.end(); ++tdii) {
      maxvals[tdii->first] = 10000.0 * normalization_factors[tdii->first];
    }

    // for each image, calculate the scores for the three features and add them
    double score;
    for (uint i=0; i<N; ++i) {
      score = 0.0;
      for (tdii=querytextdistindices.begin(); tdii!=querytextdistindices.end(); ++tdii) {
        if (distsToImages[tdii->first][i] != maxvals[tdii->first]) {
          score += exp(-distsToImages[tdii->first][i]);
        }
      }
      if (score > 0.0) {
        ResultPair r;
        r.first = score;
        r.second = i;
        result.push_back(r);
      }
    }
  }

  DBG(10) << result.size() << " results." << endl;

  return result;
}

::std::pair< ::std::vector< ::std::string >,::std::vector< ::std::string > > Retriever::getMetaFeatureInfo() {
  return database_.getMetaFeatureInfo();
}

string Retriever::dist(const uint idx, BaseDistance* dist) {
  imageComparator_.distance(idx, dist);
  ostringstream oss("");
  oss << "dist " << idx << " " << dist->name();
  return oss.str();
}

string Retriever::weight(const uint idx, const double w) {
  LinearScoring* s=dynamic_cast<LinearScoring*>(scorer_);
  ostringstream oss("");
  if (s) {
    s->weight(idx)=w;
    oss << "weight " << idx << " " << w;
  } else {
    oss << "setting weight not supported";
  }
  return oss.str();
}

string Retriever::filelist(const string filelist, string partialLoadingString) {
  DBG(10) << "Reading filelist: " << filelist << endl;
  database_.clear();
  uint nr=database_.loadFileList(filelist);
  if (nr<=0) {
    return "filelist FAILURE";
  }
  // after the filelist is processed and found to be good the partialloading is initialized
  if (partialLoadingApply_ && (partialLoadingString=="default" || partialLoadingString!="empty")) {
    this->setPartialLoading(partialLoadingString);
  } else {
    partialLoadingApply_=false;
  }
  DBG(10) << "Loading features for " << database_.size() << " images..."<< endl;
  database_.loadFeatures();
  DBG(10) << "Finished loading features." << endl;
  ostringstream oss("");
  imageComparator_=ImageComparator(database_.numberOfSuffices());
  for (uint i=0; i<database_.numberOfSuffices(); ++i) {
    imageComparator_.distance(i, new BaseDistance());
  }
  oss << "filelist " << filelist << " " << nr;
  return oss.str();
}

string Retriever::results(const uint res) {
  results_=res;
  ostringstream oss("");
  oss << "results = " << results_;
  return oss.str();
}

string Retriever::extensions(const uint ext) {
  extensions_=ext;
  ostringstream oss("");
  oss << "extensions = " << extensions_;
  return oss.str();
}

string Retriever::random(uint nr) {
  ostringstream out;
  int rndEntry;
  srand((unsigned) time(0));
  for (uint i=0; i<nr; ++i) {
    rndEntry=rand() % database_.size();
    out << database_.filename(rndEntry) << " ";
  }
  return out.str();
}

uint Retriever::results() {
  return results_;
}

BaseDistance* Retriever::dist(const uint idx) {
  return imageComparator_.distance(idx);
}

double Retriever::weight(const uint idx) const {
  LinearScoring* s=dynamic_cast<LinearScoring*>(scorer_);
  double w=0;
  if (s) {
    w=s->weight(idx);
  } else {
    ERR << "Weights not supported in this Scorer!" << endl;
  }
  return w;
}

string Retriever::filelistEntries() const {
  ostringstream oss;
  for (uint i=0; i<database_.size(); ++i) {
    oss << database_.filename(i) << " ";
  }
  return oss.str();
}

string Retriever::filelist(const uint nr) const {
  return database_.filename(nr);
}

bool Retriever::haveClasses() const {
  return database_.haveClasses();
}

int Retriever::clas(const uint &idx) const {
  const ImageContainer *t=NULL;
  t=database_[idx];
  if (t) {
    return t->clas();
  } else {
    return -1;
  }
}

int Retriever::clas(const string &filename) const {
  const ImageContainer *t=NULL;
  t=database_.getByName(filename);
  if (t) {
    return t->clas();
  } else {
    return -1;
  }
}

string Retriever::info() const {
  ostringstream oss;
  oss << "filelist " << database_.size() << " " << "results " << results_ << " " << "path " << database_.path() << " " << "extensions " << extensions_ << " " << "scoring " << scorer_->type() << " " << scorer_->settings();
  for(uint i=0;i<imageComparator_.size() and i<database_.numberOfSuffices();++i) {
    oss << " suffix "<< i << " " << database_.suffix(i) << " " <<imageComparator_.distance(i)->name();
  }
  return oss.str();
}

uint Retriever::numberOfFilelistEntries() const {
  return database_.size();
}

    void Retriever::printImageInformation(const ::std::string& imagename, ::std::ostream & os) {
      database_.printImageInformation(imagename, os);
    }

    void Retriever::printFeatureInformation(const ::std::string& imagename, const ::std::string& featureno, ::std::ostream & os) {
      istringstream is(featureno);
      uint no;
      is >> no;
      ::std::string dist_name = imageComparator_.distance(no)->name();
      if (dist_name!="textfeature") {
        database_.printFeatureInformation(imagename, featureno, os);
      } else {
        printTextFeatureInformation(imagename, featureno, os);
      }
    }

    ///print the specified text feature to the given stream. Called by printFeatureInformation
    void Retriever::printTextFeatureInformation(const ::std::string& imagename, const ::std::string& featureno, ::std::ostream & os) {

      //
      // Get the server information
      //
      istringstream is(featureno);
      uint no;
      is >> no;
      TextFeatureDistance *dist = (TextFeatureDistance*)imageComparator_.distance(no);

      ::std::string server;
      unsigned port;
      std::string language;

      dist->getServerSettings(server, port, language);

      //
      // Get the filename of the text file
      //
      ImageContainer* cont = database_.getByName(imagename);
      TextFeature *tf = (TextFeature*)cont->operator[](no);
      ::std::string textfilename = tf->value();

      //
      // Retrieve the text
      //
      ::std::string query = ":printfile "+textfilename;

      DBG(15) << "trying to contact WMIR ...";
      Socket sock(server, port);
      if (!sock.connected()) {
        DBG(15) << "not OK" << endl;
        ERR << "Could not connect! Make sure WMIR is started." << endl;
      } else {
        DBG(15) << "OK" << endl;
      }

      // Query WMIR and store the result in the rsv_table
      DBG(30) << "Sending query" << endl;
      sock << query + "\r\n";
      string text;
      DBG(30) << "getting text" << endl;
      //text = sock.getline();
      text = sock.receive();

      sock << ":bye\r\n";
      sock.close();

      //
      // Print it
      //
      DBG(10) << "Got text: " << text << endl;
      os << "suffix " << database_.suffix(no) << endl;
      os << "Textfile \": " << textfilename << "\"" << endl;
      os << text << endl;

      os << "end";
    }

    void Retriever::setFilterApply(bool apply) {
      filterApply_ = apply;
    }

    void Retriever::setFilter(vector< pair<uint,uint> > filter) {
      filter_ = filter;
    }

    void Retriever::clearFilter() {
      filter_.clear();
    }

    bool Retriever::setPartialLoadingApply(bool apply) {
      // if partial loading is to be switched on a sanity check
      // against filterApply_ has to be done
      if (apply) {
        if (filterApply_) {
          partialLoadingApply_=apply;
          return (true);
        } else {
          partialLoadingApply_=false;
          return (false);
        }
      } else {
        partialLoadingApply_=apply;
        return (false);
      }
      return (false);
    }

    void Retriever::setPartialLoading(string unparsed) {
      //sanity check against filterApply_ and filter_
      //if both aren't set or contain at least one value
      // there is no reason in performing partialLoading at all
      if (filterApply_ && filter_.size()>0) {
        uint suffSize = database_.numberOfSuffices();
        if (unparsed!="default") {
          vector<bool> dontLoad(suffSize, false);
          // parse the string
          istringstream inss(unparsed);
          while (!inss.eof() && inss.peek()!=EOF) {
            string tmp;
            // set default delimiter to our choosen one
            getline(inss, tmp, ':');
            dontLoad[atoi(tmp.c_str())]=true;
          }
          // now ensure that the bit representing the first feature selected in the filter
          // is not set
          dontLoad[filter_[0].first]=false;
          database_.setNotToLoad(dontLoad);
        } else { //default partial loading i.e. all features specified in the filter will be loaded
          vector<bool> dontLoad(suffSize, true);
          for (uint i=0; i<filter_.size(); ++i) {
            dontLoad[filter_[i].first]=false;
          }
          database_.setNotToLoad(dontLoad);
        }
      } else {
        DBG(10) << "sanity check partial loading vs. filtered retrieval failed" << endl;
        DBG(10) << "no partial loading; loading all features." << endl;
        partialLoadingApply_=false;
      }
    }

    ::std::string Retriever::getFeatures() const {
      return database_.getFeatures();
    }

    ::std::string Retriever::getPath() const {
      return database_.path();
    }
    ::std::string Retriever::getT2bPath() const {
      return database_.t2bpath();
    }

    void Retriever::loadQuery(const string &filename) {
      ImageContainer *t=new ImageContainer(filename, database_.numberOfSuffices());
      database_.loadQuery(filename, t);
      database_.addToDatabase(filename, t);

    }
