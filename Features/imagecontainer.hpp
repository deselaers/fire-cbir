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
#ifndef __imagecontainer_hpp__
#define __imagecontainer_hpp__

#include <set>
#include <vector>
#include "diag.hpp"
#include "basefeature.hpp"

typedef ::std::set< ::std::string > DescriptionSet;

/**
 * Class for an image in an image database.
 * Here we store
 *   filename
 *   features
 *   class information (optional)
 *   description information (optional)
 * and allow access to all these data.
 *
 * the class and description information is for performance evaluation
 * purposes only. If some information of this kind is to be used for
 * retrieval it has to be stored in an meta information feature and
 * the appropriate distance function has to be defined.
 */
class ImageContainer {
private:

  /// the filename of the image to which we store information
  ::std::string basename_;

  /// this is the container storing all the features to the image
  ::std::vector<BaseFeature*> features_;
  
  /// what class is this image
  uint class_;

  /// do we have additional textual information. 
  DescriptionSet description_;

  

public:
  /// default constructor. 
  ImageContainer();

  /// create an image container with enough space for numberOfFeatures
  /// features and the basename set to filename
  ImageContainer(const ::std::string& filename, uint numberOfFeatures=0);
  
  /// a copy constructor doing a deep copy if the features
  ImageContainer(const ImageContainer& src);


  /// Time to die!
  ~ImageContainer();

  /// return the basename of this image 
  const ::std::string& basename() const;

  /// return a reference to the clas of this image
  const uint& clas() const;

  /// return a reference to the clas of this image
  uint& clas();

  /// return a reference to the description of this image
  const DescriptionSet& description() const;

  /// return a reference to the description of this image
  DescriptionSet& description();

  /// return, how many features can be stored for this image
  const uint numberOfFeatures() {return features_.size();}
  
  /// return the idx-th feature
  BaseFeature*& operator[](uint idx);

  /// return the idx-th feature
  const BaseFeature* operator[](uint idx) const;
  
  /// returns itself as a vector
  ::std::vector<double> asVector();
  
  /// create an image container exacly like the one we have... but copy everythin...and I really mean copy.
  ImageContainer deepcopy();
  
  /// clears the imagecontainer
  void clear(void);

  ImageContainer operator-(const ImageContainer & img) const;

};


#endif
