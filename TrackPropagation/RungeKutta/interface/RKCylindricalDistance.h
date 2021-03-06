#ifndef RKCylindricalDistance_H
#define RKCylindricalDistance_H

#include "TrackPropagation/RungeKutta/interface/RKDistance.h"
#include "TrackPropagation/RungeKutta/interface/RKSmallVector.h"
#include "TrackPropagation/RungeKutta/interface/CylindricalState.h"

template <typename T, int N>
class RKCylindricalDistance : public RKDistance<T,N> {
public:
 
  typedef T                                   Scalar;
  typedef RKSmallVector<T,N>                  Vector;

  virtual ~RKCylindricalDistance() {}

  virtual Scalar operator()( const Vector& a, const Vector& b, const Scalar& rho) const {
      CylindricalState astate(rho,a,1.);
      CylindricalState bstate(rho,b,1.);
      return (astate.position()-bstate.position()).mag() +
	  (astate.momentum()-bstate.momentum()).mag() / bstate.momentum().mag();
  }
 
};

#endif
