#ifndef RKOneCashKarpStep_H
#define RKOneCashKarpStep_H

#include "TrackPropagation/RungeKutta/interface/RKSmallVector.h"
#include "TrackPropagation/RungeKutta/interface/RKDerivative.h"
#include "TrackPropagation/RungeKutta/interface/RKDistance.h"

#include <utility>

template <typename T, int N>
class RKOneCashKarpStep // : RKStepWithPrecision 
{
public:

  typedef T                                   Scalar;
  typedef RKSmallVector<T,N>                  Vector;

  std::pair< Vector, T> 
  operator()( Scalar startPar, const Vector& startState,
	      const RKDerivative<T,N>& deriv,
	      const RKDistance<T,N>& dist, Scalar step);
  

};

#include "TrackPropagation/RungeKutta/src/RKOneCashKarpStep.icc"

#endif
