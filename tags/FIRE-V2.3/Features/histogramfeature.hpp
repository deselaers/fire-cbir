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
#ifndef __histogramfeature_hpp__
#define __histogramfeature_hpp__

#include <vector>
#include "vectorfeature.hpp"
#include "sparsehistogramfeature.hpp"

class HistogramFeature : public VectorFeature {
protected:
  /// the counts for the bins
  ::std::vector<uint> bins_;

  /// the normalized counts for the bins (this is for faster access to
    /// histogram data) in comparison are inherited from vectorfeature
    //::std::vector<double> data_;
  
  /// how many steps in which dimension
  ::std::vector<uint> steps_;
  
  /// the min and the max value for the dimensions
  ::std::vector<double> min_, max_;

  /// the stepsize for the bins in dimensions
  ::std::vector<double> stepSize_;

  /// how many bins in total
  uint size_;
  
  /// dimensionality of input data, dim_==stepSize_.size()==steps_.size()==min_.size()==max_.size()
  uint dim_;
  
  /// how many events have been counted
  uint counter_;


public:
  /*------------------------------------------------------------
    Constructor/Destructor
    ------------------------------------------------------------*/
  /// get a histogram with size many bins
  HistogramFeature(uint size);
  
  /// get a histogram with steps[i] many steps in the i-th dimension.
  HistogramFeature(::std::vector<uint> steps);
  
  /// get a histogram with basically no space
  HistogramFeature();

  /// make a histogram containing the bins and the counter
  HistogramFeature(uint counter, const ::std::vector<uint> & bins);

  HistogramFeature(SparseHistogramFeature sh);
    
//  HistogramFeature(const HistogramFeature &src) :VectorFeature::VectorFeature(src) {
//    bins_=src.bins();
//    steps_=src.steps();
//    min_=src.min();
//    max_=src.max() ;
//    stepSize_=src.stepSize();
//    size_=(src.size());
//    dim_=(src.dim());
//    counter_=(src.counter()) ;
//  }

//  HistogramFeature &operator=(const HistogramFeature & src) {
//    (*this)=HistogramFeature(src);
//    bins_=(src.bins());
//    data_=(src.ndata());
//    steps_=src.steps();
//    min_=(src.min());
//    max_=(src.max());
//    stepSize_=(src.stepSize());
//    size_=(src.size());
//    dim_=(src.dim());
//    counter_=(src.counter());
//    return *this;
//  }

  /// destroy a histogram
  virtual ~HistogramFeature();

  virtual HistogramFeature* clone() const {
    return new HistogramFeature(*this);
  }
  /*------------------------------------------------------------
    Loading/Saving
    ------------------------------------------------------------*/
  //inherited from BaseFeature
  //virtual void load(const ::std::string &filename);
  //virtual void save(const ::std::string &filename);
  
  /// load a histogram of the old type (old FIRE versions)
  virtual void loadOld(const ::std::string &filename);

  ///inherited from BaseFeature, used in load and largefeaturefiles
  virtual bool read(::std:: istream & is);
  
  ///inherited from BaseFeature, used in loadBinary and largebinaryfeaturefiles
  virtual bool readBinary(::std::istream & is);
  
  ///inherited from BaseFeature, used in write and saving largefeaturefiles
  virtual void write(::std::ostream & os);  

  ///inherited from BaseFeature, used in writeBinary and saving largebinaryfeaturefiles
  virtual void writeBinary(::std::ostream & os);  

  /*------------------------------------------------------------
    Access to the data
    ------------------------------------------------------------*/

  ///derived from vector feature, but overwritten: How many bins in total
  virtual const uint size() const;
  
  ///derived from vector feature, but overwritten: How large is one instance in binary
  virtual const unsigned long int calcBinarySize() const;
  
  /// get the relative frequency for bin no idx
  inline virtual double operator[](const uint idx) const {return data_[idx];}
  
  /// get the relative frequency for the bin at the position pos
  virtual double operator()(const ::std::vector<uint>& pos) const;
  
  /// get the absolute count of bin idx
  virtual const uint& bin(const uint idx) const ;

  /// increase the count of bin idx
  virtual void feedbin(const uint idx);
  
  /// get the absolute count of the bin at the given position
  virtual const uint& bin(const ::std::vector<uint>& pos) const;

  /// increase the count of the bin at the given position
  virtual void feedbin(const ::std::vector<uint>& pos);

  /// increase the count of the bin, where the given points falls into
  virtual void feed(const ::std::vector<double>& point);

  /// return the internal data structure
  virtual const ::std::vector<double> &ndata() const {return data_;}

  /// return the internal data structure
  virtual const ::std::vector<uint> &bins() const {return bins_;}
  
  virtual const ::std::vector<uint> &steps() const {return steps_;}

  /// return the counter
  virtual const uint counter() const {return counter_;}
  virtual uint counter() {return counter_;}
  virtual const ::std::vector<double> &stepSize() const {return stepSize_;}
  
  const int dim() const {return dim_;}
  
  /*------------------------------------------------------------
    Helper
    ------------------------------------------------------------*/

  /// find the position vector of the bin where the given data point
  /// falls into
  ::std::vector<uint> pointToPos(const ::std::vector<double>& point) const;
  
  /// get the number of a bin, given its position vector
  uint posToBin(const ::std::vector<uint>& pos) const;

  /// initialize the stepsize, this is necessary before any "feeding"
  /// can be done.  this needs min and max to be set
  void initStepsize();

  /// set/get the maximum value to be entered into this histogram
  virtual ::std::vector<double>& max() {return max_;}
  /// get the maximum value to be entered into this histogram (const)
  virtual const ::std::vector<double>& max() const {return max_;}
  
  /// set/get the minimum value to be entered into this histogram
  virtual ::std::vector<double>& min() {return min_;}
  /// get the minimum value to be entered into this histogram
  virtual const ::std::vector<double>& min() const {return min_;}
  
};

#endif
