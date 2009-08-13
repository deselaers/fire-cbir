
/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

FIRE is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software Foundation, Inc.,
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

// NOTE: Make sure to compile this class with
// GLIBCPP_FORCE_NEW = 1
// and
// GLIBCXX_FORCE_NEW = 1
// to prevent running into memory problems !!!

#ifndef __sparsehistogramfeature_hpp__
#define __sparsehistogramfeature_hpp__

#include "basefeature.hpp"
#include "string.h"
#ifdef _SUNPRO_CC
#include <hash_map>
#else
#include <ext/hash_map>
#endif
#include <vector>

typedef std::string Position;

// comparator for strings, simply uses std::string comparison
struct PositionEquals {
  bool operator()(const Position& s1, const Position& s2) const  {
    return (s1 == s2);
  }
};

// hash function for the hash maps - simply calls builtin hash function for char*
#ifdef _SUNPRO_CC
struct PositionHash { 
  std::hash<char*> h;
  size_t operator()(const Position& p) const { return h(p.c_str()); }
};
#else
struct PositionHash { 
  __gnu_cxx::hash<char*> h;
  size_t operator()(const Position& p) const { 
    
    //    size_t conv=0;
    
    return h(p.c_str()); 

  }
};
#endif

// representation of a point
typedef std::vector<double> SHPoint;

// declaration of the hash map (storing integer or double values)
#ifdef _SUNPRO_CC
typedef std::hash_map<const Position, uint, PositionHash, PositionEquals> MapTypeInt;
typedef std::hash_map<const Position, double, PositionHash, PositionEquals>
#else
typedef __gnu_cxx::hash_map<const Position, uint, PositionHash, PositionEquals> MapTypeInt;
typedef __gnu_cxx::hash_map<const Position, double, PositionHash, PositionEquals>
#endif
MapTypeDouble;

class SparseHistogramFeature: public BaseFeature {
  
public:
  SparseHistogramFeature();
  SparseHistogramFeature(std::vector<uint> steps);
  SparseHistogramFeature(std::vector<uint> steps,const ::std::vector<double>& min, const ::std::vector<double>& max);

  SparseHistogramFeature(uint dimensions, uint steps);
  
  virtual SparseHistogramFeature* clone() const {return new SparseHistogramFeature(*this);}

  virtual ~SparseHistogramFeature();
  
  // access to the values in the bins
  double operator()(const Position& pos);
  uint binCount(const Position& pos);
  MapTypeDouble getMap() const;
  MapTypeInt getBinMap() const;

  // add items to the histogram
  void feedbin(const Position& pos);
  Position feed(const SHPoint& point);
  // must be called after the most recent call to feed/feedbin before using the () operator
  void calcRelativeFrequency();

  // get information about the size of the histogram
  virtual const uint size() const;
  // return total document length (sum of bins)
  virtual const uint length() const;
 
  // calculate the size in binary of an instance of this class
  virtual const unsigned long int calcBinarySize() const;

  const uint counter() const;

  /// set/get the maximum value to be entered into this histogram
  virtual ::std::vector<double>& max();
  /// get the maximum value to be entered into this histogram (const)
  virtual const ::std::vector<double>& max() const;
  
  /// set/get the minimum value to be entered into this histogram
  virtual ::std::vector<double>& min();
  /// get the minimum value to be entered into this histogram
  virtual const ::std::vector<double>& min() const;

  // to be called after the maximum and minimum value to be entered into the histogram have been set
  // and before any items are actually entered
  void initStepsize();

  // I/O
  virtual bool read(::std:: istream & is);
  // read binary data
  virtual bool readBinary(::std::istream &is);
  virtual void write(::std::ostream & os);
  // write histogram to given stream in binary mode
  virtual void writeBinary(::std::ostream &os);

  // smoothing of the histogram
  void expand(MapTypeDouble& expandedMap, const double smoothFactor) const;
  void clearExpandedMap(MapTypeDouble& expandedMap) const;

  // dismiss neglectable bins
  uint prunePoorBins(uint threshold);
  uint prunePoorBins(double threshold);

  // distance calculation helper methods
  MapTypeDouble getDifferenceHistogram(const MapTypeDouble& otherMap);
  void getFilledNeighborBins(const Position& pos, std::vector<Position>& neighbors);

  // create the "non-sparse" version of this histogram
  void getAllBinValues(::std::vector<double>& allBinValues) const;
  // create the "non-sparse" smoothed version of this histogram
  void getAllSmoothedBinValues(::std::vector<double>& allSmoothedBinValues, double smootheFactor) const;

  uint dimensions() { return dimensions_; }
  ::std::vector<uint> steps() { return steps_; }

private:

  // bin counter
  MapTypeInt bins_;
  // normalized bin values
  MapTypeDouble data_;

  // number of steps for each dimension
  std::vector<uint> steps_;
  // step size for each dimension
  std::vector<double> stepSize_;

  /// the min and the max value for the dimensions
  std::vector<double> min_, max_;
  
  uint dimensions_;
  uint counter_;
  uint length_;
  
  Position pointToPos(const SHPoint& point) const;

  void collectBinValues(const MapTypeDouble& theMap, ::std::vector<double>& binValues) const; 

};

#endif
