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
#ifndef __database_hpp__
#define __database_hpp__

#include <vector>
#include <string>
#include "imagecontainer.hpp"
#include "diag.hpp"
#include "featureloader.hpp"
#include "largebinaryfeaturefile.hpp" 

class Database {
private:

  /// this is the main part of this class. it does contain all ImageContainer
  ::std::vector<ImageContainer *> database_;

  /// this is for fast name based access to the data
  ::std::map< ::std::string, uint > name2IdxMap_;

  /// what suffices, that is what types of features does every image have?
  ::std::vector< ::std::string > suffixList_;
  
  /// what external files are used for binary access 
  ::std::vector<LargeBinaryFeatureFile *> suffixBinFiles_;
  
  /// stating whether a feature is not to be loaded from a largebinaryfeature
  /// at start up, or not
  /// only used in the context of filtered retrieval 
  ::std::vector<bool> binFilesNotToLoad_;
  	
  /// an object that handles the loading of features. The
  /// file suffices are mapped to the right feature type, the feature
  /// type is constructed and the file is loaded.
  FeatureLoader fl;

  /// the features are in the same directory as the images or in
  /// subdirectories of the database directory with the names of the suffices
  bool featuredirectories_;
  
  /// we have class information?
  bool classes_;
  
  /// we have a textual description for the images
  bool descriptions_;
  
  /// large feature files are used. 
  bool largefeaturefiles_;

  /// large binary feature files are used.
  bool largebinaryfeaturefiles_;
  	
  /// the basepath of the database
  ::std::string path_;

  /// the path to the type2bin file
  ::std::string t2bpath_;

  ///  a method that compares whether the two given features are
  ///  consistent. returns true if they are, false otherwise. But true
  ///  is only a "probably true"
  bool checkConsistency(const BaseFeature *ref, const BaseFeature *test) const;
  
public:


  /// constructor
  Database() : featuredirectories_(false), classes_(false), descriptions_(false), largefeaturefiles_(false), largebinaryfeaturefiles_(false), path_("") {}
  ~Database();


  /// load an image that is not in the database and the appropriate
  /// features. This method is implemented here as only the database
  /// knows which features have to be loaded. This won't use any large
  /// feature file but assumes that the features for this file are
  /// given in single files
  bool loadQuery(const ::std::string& filename, ImageContainer *result);

  /// return the idx-th ImageContainer object
  ImageContainer* operator[](uint idx) { return database_[idx]; }

  /// return the idx-th ImageContainer object
  const ImageContainer* operator[](uint idx) const { return database_[idx]; }

  /// return the idx-th LargeBinaryFeatureFile object
  LargeBinaryFeatureFile* operator()(uint idx) {return suffixBinFiles_[idx]; }
  
  /// return the idx-th LargeBinaryFeatureFile object
  const LargeBinaryFeatureFile* operator()(uint idx) const {return suffixBinFiles_[idx]; }
  
  /// search an ImageContainer given it's name. Take care, this is
  /// inefficiently implemented (linear search)
  ImageContainer* getByName(const ::std::string &filename) const;

  
  //remove the image with index i from the database
  void remove(const uint i); 

  /// return whether featuredirectories are used or not
  bool featuredirectories() const {return featuredirectories_;}

  /// clear the database. afterwards there won't be any features left
  void clear();

  /// load a filelist (but not the features)
  /// TODO: document format of filelist
  uint loadFileList(::std::string filelist);
  
  /// load the features specified by a previously loaded filelist
  void loadFeatures();
 
  /// how many images are in this database
  uint size() const { return database_.size(); }

  /// what is the filename of the idx-th image
  ::std::string filename(const uint idx) const;

  ///print all information about the specified image to the given stream
  void printImageInformation(const ::std::string& imagename, ::std::ostream & os);

  ///print all information about the specified image to the given stream
  void printFeatureInformation(const ::std::string& imagename,const ::std::string& featureno,::std::ostream & os);

  /// what is the basepath of the database. 
  const ::std::string& path() const {return path_;}
  /// what is the path of the type2bin-file
  const ::std::string& t2bpath() const {return t2bpath_;}

  /// how many suffices (that is, how many features per image) do we have?
  uint numberOfSuffices() const; 
  
  /// give the idx-th suffix
  const ::std::string& suffix(const uint idx) const { return suffixList_[idx];}
  
  /// give the relevant part of the i-th suffix
  ::std::string  relevantSuffix(const uint i) const;

  /// do we have classes?
  bool haveClasses() const;

  FeatureType featureType(const uint i) const {return fl.suffix2Type(relevantSuffix(i));}
  
  /// get a list of the available metafeature keys with example values
  ::std::pair< ::std::vector< ::std::string >,
               ::std::vector< ::std::string > > getMetaFeatureInfo() const;
  
  /// return features
  ::std::string getFeatures() const;
  /// adds image to database
  void addToDatabase(const ::std::string& filename, ImageContainer *newic);
  
  ///set the vector indicating which LargeBinaryFeatures are to be loaded
  void setNotToLoad(::std::vector<bool>& flags);
  
  /// load the feature stored in one off the lbff in suffixBinFiles_
  /// for all images in the database.
  /// return value is true if no problems where encountered while loading
  bool loadFromLBFF(const uint lbffidx);

  /// load the feature stored in one off the lbffs in suffixBinFiles_
  /// for all images specified by their number in the database
  /// retrun value is true if no problems where encountered while loading
  bool loadFromLBFF(const uint lbffidx,const ::std::vector<uint>& images);
  
  /// remove the feature information for feature with the index number idx
  /// from the database for all imageIDs in reallyLoaded
  /// this method is mainly used in the context of filtered retrieval with applied
  /// partial loading of features
  void removeFeatureInformation(uint idx, ::std::vector<uint>& reallyLoaded);

  ///return the boolean flag value for partial loading for the given feature index
  const bool binFilesNotToLoad(uint idx);
  
};

#endif
