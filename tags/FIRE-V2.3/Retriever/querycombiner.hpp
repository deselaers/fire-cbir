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
#ifndef __querycombiner_hpp__
#define __querycombiner_hpp__
#include <vector>
#include "retriever.hpp"
#include "imagecontainer.hpp"
#include "supportvectormachine.hpp"
/** this class is the prototype for relevance feedback techniques
 *  this class is abstract! 
 */
class QueryCombiner {
public:
  QueryCombiner() {}
  QueryCombiner(Retriever&) { }
  virtual ~QueryCombiner(){}
  
  virtual void query(const std::vector<ImageContainer*> posQ, 
                     const std::vector<ImageContainer*> negQ, 
                     std::vector<ResultPair>& results)=0;
  
  virtual void setParameters(const std::string&) {}
};

/** this is the default query combiner that was used in FIRE between 2003 and early 2008.
 * it is just adding up positive and negative scores to obtain a combined score
 */
class ScoreSumQueryCombiner : public QueryCombiner {
public:
  ScoreSumQueryCombiner(Retriever& r);
  virtual ~ScoreSumQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);

  virtual void setParameters(const std::string& parameters);
  
private:
  double posWeight_;
  double negWeight_;
  Retriever& retriever_;
  
};


/**
 * this is another very simple query combiner
 * we are just using the scores for each query. 
 * this method is related to the method described in NIPS 2004 by  Giorgio Giacinto and Fabio Roli.
 */
class NNQuotientQueryCombiner : public QueryCombiner {
public:
  NNQuotientQueryCombiner(Retriever& r);
  virtual ~NNQuotientQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  
private:
  Retriever& retriever_;
};

class NNScoreQueryCombiner : public QueryCombiner {
public:
  NNScoreQueryCombiner(Retriever& r);
  virtual ~NNScoreQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  
private:
  Retriever& retriever_;
};


/**
 * this is another very simple query combiner
 * we are just using the scores for each query. 
 * this is the method by Fabio Rolli and Giacinto
 */
class RelevanceScoreQueryCombiner : public QueryCombiner {
public:
  RelevanceScoreQueryCombiner(Retriever& r);
  virtual ~RelevanceScoreQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  
private:
  Retriever& retriever_;
};


class SumQuotientQueryCombiner : public QueryCombiner {
public:
  SumQuotientQueryCombiner(Retriever& r);
  virtual ~SumQuotientQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  
private:
  Retriever& retriever_;
};


class DistSumQuotientQueryCombiner : public QueryCombiner {
public:
  DistSumQuotientQueryCombiner(Retriever& r);
  virtual ~DistSumQuotientQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);
  
private:
  Retriever& retriever_;
  enum DSQMode {negDenom, allDenom};
  enum DSQMode mode_;
  


};

/** this query combiner does not use distances and scoring, but 
 * trains an SVM to decide whether images are relevant or not*/
class SVMQueryCombiner : public QueryCombiner {
public:
  SVMQueryCombiner(Retriever& r);
  virtual ~SVMQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);
private:
  Retriever& retriever_;
  int kerneltype_, degree_;
  double coef0_, gamma_, cost_;
};

class RocchioRelevanceFeedbackQueryCombiner : public QueryCombiner {
public:
  RocchioRelevanceFeedbackQueryCombiner(Retriever& r);
  virtual ~RocchioRelevanceFeedbackQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);

private:
  Retriever& retriever_;
  double alpha_, beta_, gamma_;
};

/** using weighted distance from Roberto paredes
 */
class WeightedDistanceQueryCombiner : public QueryCombiner {
public:
  WeightedDistanceQueryCombiner(Retriever& r);
  virtual ~WeightedDistanceQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);

  std::string printSigmas() const;
private:
  double WeightedL1Distance(float *v1, const std::vector<double>& v2);
  void calcWeights();
  void copiasig(float *sg, const float *sigmas);
  float wdisl1(long int i, long int j); 
  float sigmoide(float x);
  
  
  Retriever& retriever_;
  int samples, dim, classes;
  float **datos;
  int *datclas,*sel;
  float *sigmas;
  float rampa;
  float FLT_MAX, FLT_MIN;
  
  int maxIterations_;
  float stepWidth_;
  float regularisationWeight_, gradientWeight_;
  
  enum WDMode {WDRelevanceScore, WDDistSumQuotient, WDScoreSum};
  enum WDMode mode_;
};

/** class dependent weighted distance from Roberto paredes
 */
class ClassDependentWeightedDistanceQueryCombiner : public QueryCombiner {
public:
  ClassDependentWeightedDistanceQueryCombiner(Retriever& r);
  virtual ~ClassDependentWeightedDistanceQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);

  std::string printSigmas() const;
private:
  double WeightedL1Distance(float *v1, const std::vector<double>& v2, int cls);
  void calcWeights();
  void copiasig(float **sg, float **sigmas);
  float wdisl1(long int i, long int j, int sel); 
  float sigmoide(float x);
  
  Retriever& retriever_;
  int samples, dim, classes;
  float **datos;
  int *datclas,*sel;
  float **sigmas;
  float rampa;
  float FLT_MAX, FLT_MIN;
  enum WDMode {WDRelevanceScore, WDDistSumQuotient};
  enum WDMode mode_;
};

class QueryWeightingQueryCombiner: public QueryCombiner {
public:
  QueryWeightingQueryCombiner(Retriever& r);
  virtual ~QueryWeightingQueryCombiner();
  virtual void query(const std::vector<ImageContainer*> posQ, const std::vector<ImageContainer*> negQ, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);
  double distBetweenImages(const ImageContainer* q1, const ImageContainer* q2);
private:
  void calcWeights(const std::vector<ImageContainer*>& posQ, const std::vector<ImageContainer*>& negQ, std::vector<double>& posWeights, std::vector<double>& negWeights);
  float sigmoide(float x);
  
  enum WeightingMethod {WeightScores, WeightDistances};
  enum WeightingMethod mode_;
  
  Retriever& retriever_;
  float FLT_MAX, FLT_MIN;
  float rampa;
  double posWeight_, negWeight_;

};

#endif
