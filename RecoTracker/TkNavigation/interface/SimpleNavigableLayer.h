#ifndef TkNavigation_SimpleNavigableLayer_H
#define TkNavigation_SimpleNavigableLayer_H

#include "TrackingTools/DetLayers/interface/NavigableLayer.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"

#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
#include "RecoTracker/TkNavigation/interface/TkLayerLess.h"

#include <vector>

class GtfPropagator;

/** A partial implementation of the NavigableLayer
 */

class SimpleNavigableLayer : public NavigableLayer {
public:

  typedef std::vector<const DetLayer*>              DLC;
  typedef std::vector<BarrelDetLayer*>              BDLC;
  typedef std::vector<ForwardDetLayer*>             FDLC;

  SimpleNavigableLayer( const MagneticField* field,float eps,bool checkCrossingSide=true) :
    theSelfSearch(false),
    theEpsilon(eps),thePropagator(field),theCheckCrossingSide(checkCrossingSide) {}

    virtual void setInwardLinks(const BDLC&, const FDLC&, TkLayerLess sorter = TkLayerLess(outsideIn)) = 0;
  
  virtual void setAdditionalLink(DetLayer*, NavigationDirection direction=insideOut) = 0;

  void setCheckCrossingSide(bool docheck) {theCheckCrossingSide = docheck;}

  bool theSelfSearch;

  virtual std::vector< const DetLayer * > compatibleLayers (const FreeTrajectoryState &fts, 
							    PropagationDirection timeDirection,
							    int& counter)const ;

protected:

  typedef BDLC::iterator                       BDLI;
  typedef FDLC::iterator                       FDLI;
  typedef BDLC::const_iterator                 ConstBDLI;
  typedef FDLC::const_iterator                 ConstFDLI;
  typedef TrajectoryStateOnSurface             TSOS;

  float                     theEpsilon;

  mutable AnalyticalPropagator     thePropagator;
  
  bool theCheckCrossingSide;

  bool wellInside( const FreeTrajectoryState& fts, PropagationDirection dir,
		   const BarrelDetLayer* bl, DLC& result) const;

  bool wellInside( const FreeTrajectoryState& fts, PropagationDirection dir,
		   const ForwardDetLayer* bl, DLC& result) const;

  bool wellInside( const FreeTrajectoryState& fts, PropagationDirection dir,
		   ConstBDLI begin, ConstBDLI end, DLC& result) const;

  bool wellInside( const FreeTrajectoryState& fts, PropagationDirection dir,
		   const DLC& layers, DLC& result) const;

  bool wellInside( const FreeTrajectoryState& fts, PropagationDirection dir,
		   ConstFDLI begin, ConstFDLI end, DLC& result) const;

  Propagator& propagator( PropagationDirection dir) const;

  void pushResult( DLC& result, const FDLC& tmp) const;
  void pushResult( DLC& result, const BDLC& tmp) const;

  TSOS crossingState(const FreeTrajectoryState& fts,PropagationDirection dir) const;
  
};

#endif // SimpleNavigableLayer_H
