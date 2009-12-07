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
#ifndef __retriever_hpp__
#define __retriever_hpp__

#include <string>
#include <stack>
#include "diag.hpp"
#include "imagecontainer.hpp"
#include "database.hpp"
#include "imagecomparator.hpp"
#include "basescoring.hpp"
#include "linearscoring.hpp"
#include "maxentscoring.hpp"
#include "svmscoring.hpp"
#include "getscoring.hpp"
#include "distanceinteractor.hpp"


typedef ::std::pair<double,uint> ResultPair;


class Retriever;

class QueryCombiner;
class AdddingQueryCombiner;
class ReRanker;
#include "querycombiner.hpp"
#include "reranker.hpp"


/** retriever-class: this class contains the logic of the image
 retriever. Here several classes are used which together form the
 retrieval engine.
 
 - database: here the features representing the images in which we
 can search for images are managed. The database can be seen as a
 vector of ImageContainers. Each ImageContainer is a vector of
 features representing the images. This is accompanied by some
 meta-data like file names.
 
 - ImageComparator: here the feature comparison methods are
 encapsulated. The ImageComparator is a container for image
 comparison methods and it is handed a query image and database
 images and it calls the approppriate distance functions for
 comparing images. 
 
 - Scorer: this is the class combining the results from the image
 comparator into one Retrieval result. That is, given two images,
 the image comparator return a vector of distances (as many as the
 images have features representing them), these distances are
 combined into one score. This can be seen as classification into
 relevant and irrelevant images and the returning score is the
 confidence for an image to be relevant. Thus, the images with
 highest scores are returned first.
 
 - QueryCombiner: process queries at a higher level.
  The query combiner is responsible for combining the  scores of several 
  queries this is usefull if a query consists of several documents as it 
  happens e.g. in relevance feedback. Several Query combiners were implemented 
  in Valencia in January 2008 experimenting with weighting and discriminative 
  training. All of these are contained in the files querycombiner.hpp 
  and querycombiner.cpp
  QueryCombiner has access to nearly everything in Retriever and therefore is
  a place where many things could potentially be implemented
   
 
 additionally, the Retriever-class contains many helper functions
 to allow for access to data in higher layers.
 
 */
class Retriever {

private:

  /// this database object holds the database of images from which we retrieve images
  Database database_;

  /// this object knows how to compare images from the database. That
  /// is, the distance functions and the associated weightings are
  /// known here.
  ImageComparator imageComparator_;

  /// querycombiner is the object that takes care about combining different queries if several 
  /// queries are given in one step. This allows for implementation of relevance feedback mechanisms
  QueryCombiner* queryCombiner_;

  /// reranker is able to re-rank a result list, e.g. diversity based reranking with k-means clustering
  ReRanker* reRanker_;

  /// the number of results returned for a query
  uint results_;

  /// the number of images taken from a query result which is used as
  /// automatic positive relevance feedback for query expansion
  uint extensions_;

  /// the object that calculates the score from a given distance vector
  BaseScoring *scorer_;

  /// DistanceInteractor is an object that can access all distances
  /// (basically the complete distance matrix for a query) and change
  /// them in taking into account certain interactions
  DistanceInteractor interactor_;

  /// boolean indicating whether a filtered retrieval is performed or not
  bool filterApply_;

  /// boolean indication whetcher partial loading of features is performed or not
  /// this boolean is connected with filterApply_ !
  /// if filterApply_ is notset and partialLoadingApply_ is set this is an error
  bool partialLoadingApply_;

  /// a pair containing in first postion which distance should be used and
  /// in second position how many best images will be used for further
  /// retrieval. Example filter_.first = 0 and filter_.second = 1000
  /// means that the 1000 best images regarding distance 0 will be used
  ::std::vector< ::std::pair<uint,uint> > filter_;

  /**
   * given the queries, find the appropriate ImageContainers. If a
   * name is given for which no image is in the database it is tried
   * to load the image.
   * @param queryNames the names of the queries, these have to be resolved
   * @param queries here the ImageContainers are put that belong to the queryNames
   * @param newCreated here pointers to those ImageContainers are
   *     stored that are newly created, that is loaded, in this
   *     function. This is necessary to allow the desctruction of these
   *     images as soon as they are not needed anymore.
   */
  void resolveNames(const ::std::vector< ::std::string > &queryNames, ::std::vector<ImageContainer*> &queries, ::std::stack<ImageContainer*> &newCreated);
  /**
   * given vectors stillToConsider containing the best images in regard to score of the
   * preceeding filtering steps and depreciated containing the remaining images in the
   * database getBest adapts stillToConsider to a size of amount of the best images with
   * regard to the given scores. depreciated is adapted likewise
   * @param stillToConsider vector of suitable images resulting from previous filter steps
   * @param depreciated vector of nonsuitable images
   * @param scores vector of scores: size=databasesize
   * @param amount number of how many good images are wanted: is expected to be smaller than stillToConsider.size()
   */
  void getBest(::std::vector<uint> &stillToConsider, ::std::vector<uint> &depreciated, const ::std::vector<double> &scores, uint &amount);
public:

  void setCache(const std::string filename) {
    imageComparator_.setCache(filename);
    imageComparator_.openCache();
  }

  /// default constructor
  Retriever();
  
  
  //destructr
  ~Retriever();
  
  
  /// assuming database is loaded and distances are set, distances (or other things might be initialized here)
  void initialize();

  /// given a set of positive and a set of negative example image
  /// names the retrieval process is started. This function is
  /// basically a wrapper to resolve names and
  /// retrieve(vector<ImageContainer>, vector<ImageContainer>)
  void retrieve(const ::std::vector< ::std::string >& posQueries, const ::std::vector< ::std::string >& neqQueries, ::std::vector<ResultPair>& results);

  /// given a set of positive and a set of negative example
  /// ImageContainers get the results of the retrieval. For this, the
  /// right ImageComparator is necessary.
  void retrieve(const ::std::vector<ImageContainer*>& posQueries, const ::std::vector<ImageContainer*>& negQueries, ::std::vector<ResultPair>& results);

  /// start a retrieval using some meta information
  ::std::vector<ResultPair> metaretrieve(const ::std::string& query);

  /// start a retrieval using some text
  ::std::vector<ResultPair> textretrieve(const ::std::string& query);

  /// get a list of the available metafeature keys and an example value
  // for each key
  ::std::pair< ::std::vector< ::std::string >,
  ::std::vector< ::std::string > > getMetaFeatureInfo();

  /// get the distances from the given example ImageContainer q to all images in the database
  void getScores(const ImageContainer* q, ::std::vector<double> &scores);

  /// get the scores for the given distance matrix
  void getScores(::std::vector< ::std::vector<double> > &distMatrix, ::std::vector<double> &scores);

  /// get the distances for given example ImageContainer q to all wanted postive
  /// images according to current distanceID and sets the distances for all
  /// not wanted images to 1.2 times maximum distance of the wanted images
  void getDistances(const ImageContainer* q, const ::std::vector<uint>& stillToConsider, const ::std::vector<uint>& depreciated, ::std::vector< ::std::vector<double> >& distMatrix, const uint distanceID);

  /// save the distances from the given example Image (specified by
  /// the imagename) to all database images to the specified file.
  void saveDistances(::std::string imagename, ::std::string filename);

  /// set the idx-th distance in the ImageComparator used.
  ::std::string dist(const uint idx, BaseDistance* dist);

  /// return the idx-th distance in the ImageComparator used.
  BaseDistance* dist(const uint idx);

  /// set the idx-th weight in the ImageComparator used.
  ::std::string weight(const uint idx, const double w);

  /// return the idx-th weight in the ImageComparator used.
  double weight(const uint idx) const;

  ImageComparator& imageComparator() {
    return imageComparator_;
  }
  
  Database& database() {
    return database_;
  }
  
  
  
  /// return the scoring used
  BaseScoring *scorer() {
    return scorer_;
  }

  /// return a reference to the Interactor
  DistanceInteractor& interactor() {
    return interactor_;
  }

  /// set a new interactor
  void setInteractor(::std::string interactor) {
    interactor_=DistanceInteractor(interactor);
  }

  /// return the available scorings
  ::std::string availableScorings() const {
    return listScorings();
  }

  /// set a new scoring algorithm
  void setScoring(const ::std::string &scoringname);

  /// set a new query combining algorithm
  void setQueryCombiner(const std::string &queryCombinerName);

  /// set a new reRanking algorithm
  void setReranking(const std::string &rerankingName);
  
  /// how many images are in the database?
  uint numberOfFilelistEntries() const;

  /// return the nr-th filename from the filelist
  ::std::string filelist(const uint nr) const;

  /// load the filelist with the specified name
  /// the partialLoadingString is used for initializing the partial loading datastructures
  ::std::string filelist(const ::std::string filelist, ::std::string partialLoadingString="empty");

  /// reutrn a string with all names from the filelist (starting from
  /// 0 to size(filelist)-1
  ::std::string filelistEntries() const;

  /// set the number of results to be returned in all queries from now on
  ::std::string results(const uint res);

  /// return the number of results delivered after every query
  uint results();

  /// set the number of images taken as positive examples from a query
  /// result for automatic relevance feedback a.k.a. query expansion for dummies
  ::std::string extensions(const uint ext);

  /// give a string providing information on the settings of this retriever
  ::std::string info() const;

  /// give a string containing the filenames of nr images taken randomly from the database
  ::std::string random(uint nr);

  /// return whether we have class information for the database loaded currently or not
  bool haveClasses() const;

  /// return the clas of the idx-th image from the database
  int clas(const uint &idx) const;

  /// return the class of the image with the given filename from the database
  int clas(const ::std::string &filename) const;

  /// return the number of suffices in the database
  uint numberOfSuffices() const {
    return database_.numberOfSuffices();
  }

  /// return the feature type of the i-th feature of the 0-th image in
  /// the database. This should be the feature type of the i-th
  /// feature of all images.
  FeatureType featureType(uint i) const {
    return database_.featureType(i);
  }

  ///print all information about the specified image to the given stream
  void printImageInformation(const ::std::string& imagename, ::std::ostream & os);

  ///print the specified feature to the given stream
  void printFeatureInformation(const ::std::string& imagename, const ::std::string& featureno, ::std::ostream & os);

  ///print the specified text feature to the given stream. Called by printFeatureInformation
  void printTextFeatureInformation(const ::std::string& imagename, const ::std::string& featureno, ::std::ostream & os);

  ///sets the filterApply member variable
  void setFilterApply(bool apply);

  ///sets the filter to given pair vector
  void setFilter(::std::vector< ::std::pair<uint,uint> > filter);

  ///set the partialLoadingAppy member variable
  // returns a boolean indicating if partial loading is applicable or not
  bool setPartialLoadingApply(bool apply);

  /// set the flag field in the database to correct dimension and correct flags
  /// for the partial loading of features
  void setPartialLoading(::std::string unparsed);

  /// getter function for partial loading apply
  const bool partialLoadingApply() const {
    return partialLoadingApply_;
  }

  ///sets the filter pair vector to an empty pair vector
  void clearFilter();

  /// returns all used features
  ::std::string getFeatures() const;
  /// returns image-path
  ::std::string getPath() const;
  /// return type2bin-path
  ::std::string getT2bPath() const;
  /// loads new image into database
  void loadQuery(const ::std::string &filename);

};

#endif
