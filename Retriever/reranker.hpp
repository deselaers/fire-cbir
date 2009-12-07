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
#ifndef __reranker_hpp__
#define __reranker_hpp__
#include <vector>
#include "imagecontainer.hpp"
#include "em.hpp"
#include "retriever.hpp"

/** this class is the prototype for reranking techniques, returns the old ranked list without changes
 */
class ReRanker {
public:
  
  ReRanker(){}
  ReRanker(Retriever &){}
  virtual ~ReRanker(){}
  
  virtual void rerank(const std::vector<ImageContainer*>& posQueries, 
                      const std::vector<ImageContainer*>& negQueries,
                      const std::vector<ResultPair> & oldList, 
                      std::vector<ResultPair>& results);
  
  virtual void setParameters(const std::string&){}

};


/// this class clusters the best n Results and rescores them by adding an novelty bonus based on clustermemberships and seen clusters
class ClusterReranker : public ReRanker {
public:
  ClusterReranker(Retriever & retriever): retriever_(retriever){}
  virtual ~ClusterReranker(){}
  virtual void rerank(const std::vector<ImageContainer*>& posQueries,
                      const std::vector<ImageContainer*>& negQueries,
                      const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);
  
private:
  
  Retriever & retriever_;
  int nSplits_,nConsider_,nReRank_;
  double alpha_;
};

class GreedyReranker : public ReRanker {
public:
  GreedyReranker(Retriever & retriever): retriever_(retriever){}
  virtual ~GreedyReranker(){}
  virtual void rerank(const std::vector<ImageContainer*>& posQueries,
                      const std::vector<ImageContainer*>& negQueries,
                      const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);
private:
  Retriever & retriever_;
  int nConsider_,nReRank_;
  double alpha_;
};

class DPOptimisingReranker : public ReRanker {
public: 
  DPOptimisingReranker(Retriever & retriever): retriever_(retriever){}
  virtual ~DPOptimisingReranker() {}

  virtual void rerank(const std::vector<ImageContainer*>& posQueries,
                      const std::vector<ImageContainer*>& negQueries,
                      const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);

protected:
  Retriever &retriever_;
  int nConsider_,nReRank_,nOpt_;
  double alpha_;
  bool scaleSimilarity_, noveltyMinApprox_, noveltyOnlyLast_;
};


class DPOpt2 : public DPOptimisingReranker {
public: 
  DPOpt2(Retriever & retriever): DPOptimisingReranker(retriever) {}
  virtual ~DPOpt2() {}
  
  virtual void rerank(const std::vector<ImageContainer*>& posQueries,
                      const std::vector<ImageContainer*>& negQueries,
                      const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results);
};

class VisAVisReranker : public ReRanker {
public:
  VisAVisReranker(Retriever & retriever): retriever_(retriever){}
  virtual ~VisAVisReranker(){}
  virtual void rerank(const std::vector<ImageContainer*>& posQueries,
                      const std::vector<ImageContainer*>& negQueries,
                      const std::vector<ResultPair> & oldList, std::vector<ResultPair>& results);
  virtual void setParameters(const std::string& parameters);
private:
  Retriever & retriever_;
  int nConsider_,nReRank_;
  double alpha_;
};



#endif
