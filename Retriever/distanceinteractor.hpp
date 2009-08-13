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

/** this file realizes functionality for distances to interact. 
    Several classes are defined here:

    Interaction: abstract base class which can be used to define interactions. 
    
    Each Interaction gets the distances for an image and can do to
    them whatever it wants. 

    derived from Interaction is e.g. MinInteraction which looks at a
    predefined set of distances in the vector and replaces all of them
    by the minimum of their values.

    derived from Interaciton is e.g. MaxInteraction which is similar
    to MinInteraction but takes the Maximum

    another one is MulInteraction, which replaces a set of predefined
    distances by their product.
    
    These classes are the actions that are called from the
    DistanceInteractor-class. The distance interactor class receives a
    distance matrix and applies a set of predefined interactions to
    its lines.

    In principal arbitrary many interactions can be defined to be
    applied, although I doubt that many different interactions will be
    used at the same time.

    The interactions to be applied are specified from the command
    line. An example is stated in the help of fire.
*/

#ifndef __distanceinteractor_hpp__
#define __distanceinteractor_hpp__

#include <string>
#include <vector>
#include <limits>


/** abstract base class for an interaction. given a distance vector it
    can modifiy its distances at will */
class Interaction {
public:
  Interaction() {}

  /** the vector that is given here, can e.g. specify the positions of
      the distances that are to be considered and changed. */
  Interaction(const ::std::vector<uint>&) {}

  virtual ~Interaction() {}
  
  /** this is the function that actually should analyse and modify
      distances */
  virtual void apply(::std::vector<double> & d)=0;
};


/** simple implementation of Interaction. Take the distances at the
    specified positions and replace them by their minimum*/
class MinInteraction : public Interaction { 
private: 
  ::std::vector<uint> positions_; 

public:
  MinInteraction(const ::std::vector<uint>& pos) :positions_(pos) {
  }
  
  virtual ~MinInteraction() {}

  
  virtual void apply(::std::vector<double>&d) {
    double min=::std::numeric_limits<double>::max();
    
    for(::std::vector<uint>::const_iterator i=positions_.begin();i!=positions_.end();++i) {
      if(min>d[*i]) {
        min=d[*i];
      }
    }
    
    for(::std::vector<uint>::const_iterator i=positions_.begin();i!=positions_.end();++i) {
      d[*i]=min;
    }
  }
};

/** simple implementation of Interaction. Take the distances at the
    specified positions and replace them by their maximum*/
class MaxInteraction : public Interaction {
private:
  ::std::vector<uint> positions_;

public:
  MaxInteraction(const ::std::vector<uint>& pos) :positions_(pos) {
  }
  
  virtual ~MaxInteraction() {}

  
  virtual void apply(::std::vector<double>&d) {
    double max=-::std::numeric_limits<double>::max();
    
    for(::std::vector<uint>::const_iterator i=positions_.begin();i!=positions_.end();++i) {
      if(max<d[*i]) {
        max=d[*i];
      }
    }
    
    for(::std::vector<uint>::const_iterator i=positions_.begin();i!=positions_.end();++i) {
      d[*i]=max;
    }
  }
};

/** simple implementation of Interaction. Take the distances at the
    specified positions and replace them by their product*/
class MulInteraction : public Interaction {
private:
  ::std::vector<uint> positions_;

public:
  MulInteraction(const ::std::vector<uint>& pos) :positions_(pos) {
  }
  
  virtual ~MulInteraction() {}

  
  virtual void apply(::std::vector<double>&d) {
    double mul=1;
    
    for(::std::vector<uint>::const_iterator i=positions_.begin();i!=positions_.end();++i) {
      mul*=d[*i];
    }
    
    for(::std::vector<uint>::const_iterator i=positions_.begin();i!=positions_.end();++i) {
      d[*i]=mul;
    }
  }
};


/**
   DistanceInteractor-class. The distance interactor class receives a
   distance matrix and applies a set of predefined interactions to
   its lines.
   
   In principal arbitrary many interactions can be defined to be
   applied, although I doubt that many different interactions will be
   used at the same time.
   
   The interactions to be applied are specified from the command
   line. An example is stated in the help of fire.
*/
class DistanceInteractor {
private:
  
  /// store the interactions to be applied
  ::std::vector<Interaction*> interactions_;
  
  /// given the name of an interaction extract the positions specified
  /// in the parentheses. e.g. min(2,3,4) returns the vector [2,3,4]
  ::std::vector<uint> getPositions(const ::std::string& act) {
    uint spos=act.find("(");
    uint epos=act.find(")");
    uint pos=spos+1;
    ::std::vector<uint>positions;
    while(pos<epos) {
      ::std::string s="";
      while(act[pos]!=',' and pos<epos) {
        s+=act[pos];
        ++pos;
      }
      DBG(10) << VAR(s) << ::std::endl;
      ++pos;
      positions.push_back(atoi(s.c_str()));
    }
    return positions;
  }
  
  /// given a string specifying an interaction, return the
  /// corresponding interaction
  Interaction* createInteraction(const ::std::string& act) {
    ::std::string prefix;
    prefix.assign(act,0,3);
    DBG(10) << VAR(prefix) << ::std::endl;
    
    ::std::vector<uint> positions=getPositions(act);
    
    Interaction* result=NULL;

    if(prefix=="min") {
      result=new MinInteraction(positions);
    } else if (prefix=="max") {
      result=new MaxInteraction(positions);
    } else if (prefix=="mul") {
      result=new MulInteraction(positions);
    } else {
      DBG(10) << "unknown prefix: "<< prefix << ::std::endl;
    }
    return result;
  }
  
public:
  /// default constructor: necessary for instantiation in default
  /// constructor of retriever
  DistanceInteractor() : interactions_() {}

//  ~DistanceInteractor() {
//    for(uint i=0;i<interactions_.size();++i) {
//      delete interactions_[i];
//    }
//    interactions_.clear();
//  }

  /// constructor that parses an interaction string and adopts the
  /// specified settings.
  DistanceInteractor(const ::std::string actions) {
    uint pos=0;
    ::std::string action;
    Interaction* interact;
    DBG(10) << VAR(actions) << ::std::endl;
    while(pos<actions.size()) {
      uint spos=actions.find("ACTION=",pos);
      uint epos=actions.find(":",pos+1);
      DBG(25) << VAR(spos) << VAR(epos) << VAR(pos) << ::std::endl;
      if (spos<actions.size()) {
        if(epos >= actions.size()) {epos=actions.size();}
        action.assign(actions,spos+7,epos-spos);
        DBG(10) << VAR(action) << ::std::endl;
        interact=createInteraction(action);
        interactions_.push_back(interact);
        pos=epos;
      }
    }
  }
  

  /// apply all defined actions to a given distance matrix.
  void apply(::std::vector< ::std::vector<double> >& d) {
    uint N=d.size();
    //uint M=d[0].size();
    
    for(uint i=0;i<interactions_.size();++i) {
      DBG(10) << "Applying interaction " << i << ::std::endl;
      Interaction* I=interactions_[i];
      for(uint n=0;n<N;++n) {
        
        I->apply(d[n]);
      }
    }
  }

};
#endif
