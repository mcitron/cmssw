#ifndef DATAFORMATS_TRAJECTORYSEED_TRAJECTORYSEED_h
#define DATAFORMATS_TRAJECTORYSEED_TRAJECTORYSEED_h

#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include <utility>
#include <algorithm>

/**
   TrajectorySeed contains
   - a TSOS
   - a vector of RecHits (with Own_vector to store polimorphic)
   - a propagation direction
**/
class TrajectorySeed {
 public:
  typedef edm::OwnVector<TrackingRecHit> recHitContainer;
  typedef recHitContainer::const_iterator const_iterator;
  typedef std::pair<const_iterator,const_iterator> range;
  
  
  TrajectorySeed(){}
  virtual ~TrajectorySeed(){}

  TrajectorySeed(PTrajectoryStateOnDet const & ptsos, 
		 recHitContainer const & rh, 
		 PropagationDirection  dir) : 
    hits_(rh),  tsos_(ptsos), dir_(dir) {}
  
  void swap(PTrajectoryStateOnDet & ptsos, 
		 recHitContainer & rh, 
	    PropagationDirection & dir) { 
    hits_.swap(rh);
    std::swap(tsos_,ptsos);
    std::swap(dir_,dir);
  }

  void swap(TrajectorySeed & rh) { 
    hits_.swap(rh.hits_);
    std::swap(tsos_,rh.tsos_);
    std::swap(dir_,rh.dir_);
  }

  TrajectorySeed(TrajectorySeed const & o) :
    hits_(o.hits_),  tsos_(o.tsos_), dir_(o.dir_) {}

  TrajectorySeed & operator=(TrajectorySeed const & o) {
    TrajectorySeed tmp(o); swap(tmp);   
    return *this;
  }


#if defined( __GXX_EXPERIMENTAL_CXX0X__)
  TrajectorySeed(TrajectorySeed && o) :
    hits_(std::move(o.hits_)),  tsos_(std::move(o.tsos_)), dir_(std::move(o.dir_)) {}

  TrajectorySeed & operator=(TrajectorySeed && o) {
    swap(o);   
    return *this;
  }

#endif

  range recHits() const {
    return std::make_pair(hits_.begin(), hits_.end());
  }
  unsigned int nHits() const {return hits_.size();}
  PropagationDirection direction() const {return  dir_;}
  PTrajectoryStateOnDet const & startingState() const {return tsos_;}

  virtual TrajectorySeed * clone() const {return new TrajectorySeed( * this); }

 
 private:
  edm::OwnVector<TrackingRecHit> hits_;
  PTrajectoryStateOnDet tsos_;
  PropagationDirection dir_;
};

inline void swap(TrajectorySeed & rh, TrajectorySeed & lh) { 
  rh.swap(lh);
}

typedef TrajectorySeed BasicTrajectorySeed;

#endif
