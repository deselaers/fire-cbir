#include "reranker.hpp"
#include "em.hpp"
#include "distancemaker.hpp"
#include "stringparser.hpp"
#include "basescoring.hpp"
#include "imagecomparator.hpp"

using namespace std;

void ReRanker::rerank(const std::vector<ImageContainer*>& posQueries, const std::vector<ImageContainer*>& negQueries, const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results){
  
  results=oldList;

}



bool comp(ResultPair i, ResultPair j){return i.first>j.first;}

void ClusterReranker::rerank(const std::vector<ImageContainer*>& posQueries, const std::vector<ImageContainer*>& negQueries, const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results){

  EM em;
  em.saveBeforeSplits()=false;
  em.maxSplit()=nSplits_;
  em.poolMode()=clusterPooling;
  em.splitMode()=allSplit;
  //em.stopWithNClusters()=10;
  em.disturbMode()=meanDisturb;
  em.iterationsBetweenSplits()=10;
  em.minObservationsPerCluster()=1;
  em.epsilon()=1e-5;
  em.dontSplitBelow()=1;
  DistanceMaker dm;
  em.setDist(dm.makeDistance("euclidean"));

  vector<ResultPair> imagesToConsider(nConsider_);
  results=oldList;
  sort(results.begin(),results.end(),comp);

  DoubleVectorVector trainingdata(nConsider_);
  for (int i=0;i<nConsider_;++i){
    trainingdata[i]=new DoubleVector((*retriever_.database()[results[i].second]).asVector());
  }
  
  ResultVector clusterinformation;
  
  em.run(trainingdata,clusterinformation);

  int nClusters=em.numberOfClusters();

  // now adjust scores based on cluster memberships
  IntVector clusterCounts(nClusters,0);
  int totalClusterCount=0;
  
  for (int i=0;i<nConsider_;++i){
    double noveltyBonus=alpha_;
    if (i!=0){
      noveltyBonus*=(1.0-double(clusterCounts[clusterinformation[i]])/(totalClusterCount*nClusters));
    }
    results[i].first+=noveltyBonus;
    clusterCounts[clusterinformation[i]]++;    
    totalClusterCount++;
  }
}
 


void ClusterReranker::setParameters(const std::string& parameters){
    
  nConsider_=getIntAfter(parameters,"NCONSIDER=",100);
  alpha_=getDoubleAfter(parameters,"ALPHA=",1.0);
  nSplits_=getIntAfter(parameters,"SPLITS=",5);
}

void GreedyReranker::rerank(const std::vector<ImageContainer*>& posQueries, const std::vector<ImageContainer*>& negQueries, const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results){
  results=oldList;
  double baseScore=0.0;
  if (nConsider_<int(oldList.size())){
    baseScore=oldList[nConsider_+1].first;
  }
  
  sort(results.rbegin(),results.rend());

  vector<ResultPair> bestImages(nConsider_);
  // first image is fixed
  bestImages[0]=results[0];
  bestImages[0].first=baseScore+1.0;
  
  std::vector< std::vector<double> > dissimilarityMatrix(nConsider_, std::vector<double>(nConsider_,0.0));
  
  BaseScoring* scorer=retriever_.scorer();
  ImageComparator &comparator=retriever_.imageComparator();
  
  /// CALCULATE similarity scores for the candidate images
  std::vector<double> similarityScores(nConsider_,0.0);
  double maxScore=0.0;
  for(int i=0;i<nConsider_;++i) {
    similarityScores[i]=-log(results[i].first);
    if (similarityScores[i]>maxScore) maxScore=similarityScores[i];
    DBG(20) << VAR(results[i].first) << " " << VAR(results[i].second) << " "<< VAR(similarityScores[i])<<endl;
  }
  DBG(20)<<VAR(maxScore)<<endl;
  
  for(int i=0;i<nConsider_;++i) { similarityScores[i]=1.0-similarityScores[i]/maxScore; }
  // NOW similarity scores are normalised, best one is 1, worst one is smaller than that

  /// CALCULATE pairwise dissimilarities between the images and normalise these between 1 (highest distance) and 0 (comparison of image with itself)
  double maxS=-std::numeric_limits<double>::max();
  for(int i=0;i<nConsider_;++i) {
    ImageContainer* I=retriever_.database()[results[i].second];
    for(int j=i+1;j<nConsider_;++j) {
      ImageContainer* J=retriever_.database()[results[j].second];
      double tmp=-log(scorer->getScore(comparator.compare(I,J)));
      if (tmp>maxS) maxS=tmp;
      dissimilarityMatrix[i][j]=tmp;
    }
  }
  for(int i=0;i<nConsider_;++i) {
    for(int j=i+1;j<nConsider_;++j) {
      dissimilarityMatrix[i][j]=(dissimilarityMatrix[i][j])/(maxS);
      dissimilarityMatrix[j][i]=dissimilarityMatrix[i][j];
    }
  }

  DBG(10)<<"New Order: 0";
  //find best nReRank images
  for (int i=1;i<nConsider_;++i){
    //search for i-th next best image in the @nConsider_ most relevant images
    DBG(20)<<i<<std::endl;

    bestImages[i].first=-1;
    double bestScore=-1;
    int idx=-1;
    //0-th image is ommited because it is already fixed at pos 0
    for (uint k=1;k<nConsider_;++k){
      uint itsIndex=results[k].second;
      bool alreadyFound=false;
      //look whether the image currently considered is already in the sequence up to the current position
      for (int j=0;j<i ;++j){
        if (itsIndex==bestImages[j].second){
          alreadyFound=true;
        }
      }

      if (!alreadyFound){
        double fullScore=alpha_*similarityScores[k];
        double noveltyScore=0.0;
        for (int j=0;j<i ;++j){
          noveltyScore+=dissimilarityMatrix[k][j]/i;        
        }
        DBG(20)<< similarityScores[k] <<" " << noveltyScore << endl;        
        fullScore+=(1.0-alpha_)*noveltyScore;
        if (fullScore > bestScore){
          bestScore=fullScore;
          idx=k;
        }
      }
    }
    bestImages[i].first=baseScore+1.0/double(i+1);
    bestImages[i].second=results[idx].second;
    BLINK(10)<<" "<<idx;
  }
  DBG(10)<<endl;
  for (int i=0;i<nConsider_;++i){
    results[i]=bestImages[i];
  }
}

void GreedyReranker::setParameters(const std::string& parameters){
  //nReRank_=getIntAfter(parameters,"NRERANK=",100);
  nConsider_=getIntAfter(parameters,"NCONSIDER=",100);
  alpha_=getDoubleAfter(parameters,"ALPHA=",0.8);
}


void DPOptimisingReranker::rerank(const std::vector<ImageContainer*>& posQueries, const std::vector<ImageContainer*>& negQueries, 
                                  const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results) {
  
  std::vector< std::vector<double> > dissimilarityMatrix(nConsider_, std::vector<double>(nConsider_,0.0));
  
  BaseScoring* scorer=retriever_.scorer();
  ImageComparator &comparator=retriever_.imageComparator();
  
  results=oldList;
  sort(results.rbegin(),results.rend());
  
  double basescore=0.0;
  if (nConsider_<results.size()){    basescore=results[nConsider_].first;  }
  

  /// CALCULATE similarity scores for the candidate images
  std::vector<double> similarityScores(nConsider_,0.0);
  double maxScore=0.0;
  for(int i=0;i<nConsider_;++i) {
    similarityScores[i]=-log(results[i].first);
    if (similarityScores[i]>maxScore) maxScore=similarityScores[i];
    DBG(20) << VAR(results[i].first) << " " << VAR(results[i].second) << " "<< VAR(similarityScores[i])<<endl;
  }
  DBG(20)<<VAR(maxScore)<<endl;
  
  for(int i=0;i<nConsider_;++i) { similarityScores[i]=1.0-similarityScores[i]/maxScore; }
  // NOW similarity scores are normalised, best one is 1, worst one is smaller than that

  /// CALCULATE pairwise dissimilarities between the images and normalise these between 1 (highest distance) and 0 (comparison of image with itself)
  double maxS=-std::numeric_limits<double>::max();
  for(int i=0;i<nConsider_;++i) {
    ImageContainer* I=retriever_.database()[results[i].second];
    for(int j=i+1;j<nConsider_;++j) {
      ImageContainer* J=retriever_.database()[results[j].second];
      double tmp=-log(scorer->getScore(comparator.compare(I,J)));
      if (tmp>maxS) maxS=tmp;
      dissimilarityMatrix[i][j]=tmp;
    }
  }
  for(int i=0;i<nConsider_;++i) {
    for(int j=i+1;j<nConsider_;++j) {
      dissimilarityMatrix[i][j]=(dissimilarityMatrix[i][j])/(maxS);
      dissimilarityMatrix[j][i]=dissimilarityMatrix[i][j];
    }
  }

  /// Debug output
  DBG(20) << "Similarity scores:"  ;   for(int i=0;i<nConsider_;++i) { BLINK(20) << " " << similarityScores[i] ;} BLINK(20) << endl;
  DBG(20) << "Dissimilarity Matrix:" << endl;
  for(int i=0;i<nConsider_;++i) {
    for(int j=0;j<nConsider_;++j) {
      BLINK(20) << dissimilarityMatrix[i][j] << " ";
    } 
    BLINK(20) << endl;
  }

 
  // all indices are [t][u]([predecessors])
  std::vector< std::vector<double> > H(nConsider_+1,std::vector<double>(nConsider_+1,0.0));
  std::vector< std::vector< std::vector< int > > > B2(nConsider_+1, std::vector< std::vector<int> >(nConsider_+1, std::vector<int>(0)));
  
  // init
  H[0][1]=alpha_*similarityScores[0];
  //  B2[0][1].push_back(0);
  //  for(uint i=0;i<=nConsider_;++i) {
  //   B2[0][i].push_back(0);
  //}
  

  // t=an welcher stelle in der resultliste
  // v= welcher vorgaenger wird gerade betrachtet
  // u= welches bild wird gerade als neuer kandidat betrachtet
  for(uint t=1;t<=nConsider_;++t) {
    for(uint u=1;u<=nConsider_;++u) {
      
      double bestHTUV=-std::numeric_limits<double>::max();
      int bestV=-1;
      
      double htuvSimilarity=similarityScores[u-1];

      for(uint v=1;v<=nConsider_;++v) {
        if(u==v) continue;

        
        // calculating the local score
        double htuvNovelty=0.0;
        int nFound=B2[t-1][v].size();
        
        if(noveltyMinApprox_) {
          DBG(20) << VAR(t) << " " << VAR(u) << " " << VAR(v) << " " << VAR(nFound) << " B2" ;
          DBGI(20,{for(uint i=0;i<nFound;++i) {          BLINK(10) << " " << B2[t-1][v][i] ;} BLINK(10) << endl;})
          double minNovelty=std::numeric_limits<double>::max();
          for(uint i=0;i<nFound;++i) {
            if( B2[t-1][v][i] == u-1) { // if this image is somewhere in the sequence of retrieved images before, set similarity and novelty to 0.0, we do NOT want duplicates in the results
              htuvNovelty=0.0;
              htuvSimilarity=0.0;
              break;
            }
            // dissimilarity matrix is now symmetric
            if(dissimilarityMatrix[u-1][B2[t-1][v][i]] < minNovelty) {
              minNovelty=dissimilarityMatrix[u-1][B2[t-1][v][i]];
            }
          }
          htuvNovelty=nFound*minNovelty;
        } else if(noveltyOnlyLast_) {
          htuvNovelty=dissimilarityMatrix[u-1][B2[t-1][v][nFound-1]];
          for(uint i=0;i<nFound;++i) {
            if( B2[t-1][v][i] == u-1) { // if this image is somewhere in the sequence of retrieved images before, set similarity and novelty to 0.0, we do NOT want duplicates in the results
              htuvNovelty=0.0;
              htuvSimilarity=0.0;
              break;
            }
          }

        } else {
          DBG(20) << VAR(t) << " " << VAR(u) << " " << VAR(v) << " " << VAR(nFound) << " B2" ;
          DBGI(20,{for(uint i=0;i<nFound;++i) {          BLINK(10) << " " << B2[t-1][v][i] ;} BLINK(10) << endl;})
            
            for(uint i=0;i<nFound;++i) {
              if( B2[t-1][v][i] == u-1) { // if this image is somewhere in the sequence of retrieved images before, set similarity and novelty to 0.0, we do NOT want duplicates in the results
                htuvNovelty=0.0;
                htuvSimilarity=0.0;
                break;
              }
              // dissimilarity matrix is now symmetric
              htuvNovelty+=dissimilarityMatrix[u-1][B2[t-1][v][i]];
            }
          
          // unclear why this poses problems!
          //	DBG(10) << VAR(htuvNovelty) << " " << VAR(nFound) <<  endl;
          //if(nFound>0 and htuvNovelty>0.0) {htuvNovelty/=double(nFound);}
          //	DBG(10) << VAR(htuvNovelty) << endl;
          //	double norm=1.0;
          //	if(nFound>0) norm/=nFound;
          //	DBG(10) << VAR(norm) << endl;
          //	htuvNovelty*=norm;
          if(scaleSimilarity_ and nFound>0) {htuvSimilarity*=nFound;}
        } // else 


        double htuv=(1.0-alpha_)*htuvNovelty + alpha_*htuvSimilarity;
        DBG(20) << VAR(htuv) << endl;
        if(H[t-1][v]+htuv > bestHTUV) {
          bestHTUV=H[t-1][v]+htuv;
          bestV=v;
        }
      }
      H[t][u]=bestHTUV;
      B2[t][u]=B2[t-1][bestV]; B2[t][u].push_back(bestV-1);

      DBG(20) << VAR(u) << " " << VAR(t) << " " VAR(bestHTUV) << " " << VAR(bestV) << " " << VAR(B2[t][u].size()) <<  " " << VAR(B2[t-1][bestV].size()) << " B2: ";
      for(uint i=0;i<B2[t-1][bestV].size();++i) {BLINK(20) << " " <<  B2[t-1][bestV][i];} BLINK(20) << endl;
    }
  }



  // "tracing back": does not happen, sequence is in B2
  int t=nConsider_;
  int u=max_element(H[t].begin(), H[t].end())-H[t].begin();
  u-=1;

  DBG(10) << "B2 ";
  for(uint i=0;i<B2[t][u].size();++i) {    BLINK(10) << B2[t][u][i] << " " ;} BLINK(10) << endl;
  
  std::vector<ResultPair> copyResults=results;
  
  std::vector<int> occurrenceCounter(nConsider_,0);
  for(uint i=0;i<B2[t][u].size();++i) {
    results[i].first=basescore+1.0/double(i+1);
    results[i].second=copyResults[B2[t][u][i]].second;
    occurrenceCounter[B2[t][u][i]]++;
    DBG(20) << VAR(i) << " "<< VAR(B2[t][u][i]) << " "  << VAR(results[i].first) << " " << VAR(results[i].second) << endl;
  }
  
  int count=0;
  for(uint i=0;i<nConsider_;++i){
    if (occurrenceCounter[i]==0){
      results[i].first=basescore+1.0/double(nOpt_+count+1);
      results[i].second=copyResults[i].second;
      count++;
    }
  }

  // traceback not necessary, B2 does have all information

//  vector<int> traceback;
//  
//  while(t>0) {
//    DBG(20) << VAR(u) << " " << VAR(t) << " " << VAR(H[t][u])  << endl;
//    traceback.push_back(u);
//    results[t-1].second=oldList[u-1].second;
//    u=B[t][u];
//    results[t-1].first=basescore+1.0/double(t+0.001); // really terrible way of generating a valid score
//    t-=1;
//  }
//
//  DBG(10) << "traceback ";
//  for(int i=traceback.size()-1;i>=0;--i) {    BLINK(10) << traceback[i] << " " ;} BLINK(10) << endl;
//

}



void DPOptimisingReranker::setParameters(const std::string& parameters){
  //nReRank_=getIntAfter(parameters,"NRERANK=",100);
  nConsider_=getIntAfter(parameters,"NCONSIDER=",100);
  alpha_=getDoubleAfter(parameters,"ALPHA=",0.8);
  scaleSimilarity_=getBooleanString(parameters,"SCALESIMILARITY");
  noveltyMinApprox_=getBooleanString(parameters,"NOVELTYMINAPPROX");
  noveltyOnlyLast_=getBooleanString(parameters,"NOVELTYONLYLAST");
  nOpt_=getIntAfter(parameters,"NOPT=",20);

  DBG(10) << VAR(alpha_) << " "  
          << VAR(nConsider_) << " " 
          << VAR(nOpt_) << " " 
          << VAR(scaleSimilarity_) << " " 
          << VAR(noveltyMinApprox_) << " " 
          << VAR(noveltyOnlyLast_) << endl;
}


void DPOpt2::rerank(const std::vector<ImageContainer*>& posQueries, const std::vector<ImageContainer*>& negQueries, 
                    const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results) {
  
  std::vector< std::vector<double> > dissimilarityMatrix(nConsider_, std::vector<double>(nConsider_,0.0));
  
  BaseScoring* scorer=retriever_.scorer();
  ImageComparator &comparator=retriever_.imageComparator();
  
  results=oldList;
  sort(results.rbegin(),results.rend());
  
  double basescore=0.0;
  if (nConsider_<results.size()){    basescore=results[nConsider_].first;  }
  

  /// CALCULATE similarity scores for the candidate images
  std::vector<double> similarityScores(nConsider_,0.0);
  double maxScore=0.0;
  for(int i=0;i<nConsider_;++i) {
    similarityScores[i]=-log(results[i].first);
    if (similarityScores[i]>maxScore) maxScore=similarityScores[i];
    DBG(20) << VAR(results[i].first) << " " << VAR(results[i].second) << " "<< VAR(similarityScores[i])<<endl;
  }
  DBG(20)<<VAR(maxScore)<<endl;
  
  for(int i=0;i<nConsider_;++i) { similarityScores[i]=1.0-similarityScores[i]/maxScore; }
  // NOW similarity scores are normalised, best one is 1, worst one is smaller than that
  
  /// CALCULATE pairwise dissimilarities between the images and normalise these between 1 (highest distance) and 0 (comparison of image with itself)
  double maxS=-std::numeric_limits<double>::max();
  for(int i=0;i<nConsider_;++i) {
    ImageContainer* I=retriever_.database()[results[i].second];
    for(int j=i+1;j<nConsider_;++j) {
      ImageContainer* J=retriever_.database()[results[j].second];
      double tmp=-log(scorer->getScore(comparator.compare(I,J)));
      if (tmp>maxS) maxS=tmp;
      dissimilarityMatrix[i][j]=tmp;
    }
  }
  for(int i=0;i<nConsider_;++i) {
    for(int j=i+1;j<nConsider_;++j) {
      dissimilarityMatrix[i][j]=(dissimilarityMatrix[i][j])/(maxS);
      dissimilarityMatrix[j][i]=dissimilarityMatrix[i][j];
    }
  }

  /// Debug output
  DBG(20) << "Similarity scores:"  ;   for(int i=0;i<nConsider_;++i) { BLINK(20) << " " << similarityScores[i] ;} BLINK(20) << endl;
  DBG(20) << "Dissimilarity Matrix:" << endl;
  for(int i=0;i<nConsider_;++i) {
    for(int j=0;j<nConsider_;++j) {
      BLINK(20) << dissimilarityMatrix[i][j] << " ";
    } 
    BLINK(20) << endl;
  }

 
  // all indices are [t][u]([predecessors])
  std::vector< std::vector<double> > H(nConsider_+1,std::vector<double>(nConsider_+1,0.0));
  std::vector< std::vector< int > > B(nConsider_+1, std::vector< int >(nConsider_+1));
  
  // init
  H[0][1]=alpha_*similarityScores[0];
  //  B2[0][1].push_back(0);
  //  for(uint i=0;i<=nConsider_;++i) {
  //   B2[0][i].push_back(0);
  //}
  
  
  // t=an welcher stelle in der resultliste
  // v= welcher vorgaenger wird gerade betrachtet
  // u= welches bild wird gerade als neuer kandidat betrachtet

  for(uint t=1;t<20;++t) {
    for(uint u=2;u<=nConsider_;++u) {
      
      double bestHTUV=-std::numeric_limits<double>::max();
      int bestV=-1;
      
      double htuvSimilarity=similarityScores[u-1];
      
      for(uint v=1;v<u;++v) {
        // calculating the local score
        double htuvNovelty=0.0;
        int nFound=t;
        
        htuvNovelty=dissimilarityMatrix[u-1][B[t-1][v]-1];

        double htuv=(1.0-alpha_)*htuvNovelty + alpha_*htuvSimilarity;
        if(H[t-1][v]+htuv > bestHTUV) {
          bestHTUV=H[t-1][v]+htuv;
          bestV=v;
        }
        H[t][u]=bestHTUV;
        B[t][u]=bestV;
        DBG(20)  << VAR(t) << " " << VAR(u) << " " << VAR(bestV) << endl;
      } //for v
    } // for u
  } // for t

  // traceback
  int t=19;
  int u=max_element(H[t].begin(), H[t].end())-H[t].begin();
  u-=1;
  
  std::vector<ResultPair> copyResults=results;

  vector<int> traceback;
  while(t>=0) {
    traceback.insert(traceback.begin(),u);
    u=B[t][u];
    --t;
  }
  // now the traceback array contains the top 20 in inverse ordering

  DBG(10) << "traceback ";
  for(int i=0;i<traceback.size();++i) { 
    BLINK(10) << traceback[i] << " " ;
    results[i].first=basescore+1.0/(double(i+1));
    results[i].second=copyResults[traceback[i]-1].second;
  } 
  BLINK(10) << endl;

  // for comparison
//  for(int i=0;i<20;++i) {
//    DBG(10) << VAR(i) << VAR(results[i].second)  << " "  << VAR(copyResults[i].second) << endl;
//  }
//  
  
} // DPOpt2::rerank
