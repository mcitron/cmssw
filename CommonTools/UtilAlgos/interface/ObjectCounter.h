#ifndef UtilAlgos_ObjectCounter_h
#define UtilAlgos_ObjectCounter_h
/** \class ObjectCounter
 *
 * Counts the number of objects in a collection and prints a
 * summary report at the end of a job.
 * 
 * Template parameters:
 * - C : collection type
 *
 * \author Luca Lista, INFN
 *
 * \version $Revision: 1.9 $
 *
 * $Id: ObjectCounter.h,v 1.9 2007/12/31 18:43:14 ratnik Exp $
 *
 */
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/TypeID.h"
#include <iostream>
#include <cmath>

template<typename C>
class ObjectCounter : public edm::EDAnalyzer {
public:
  /// constructor from parameter set
  explicit ObjectCounter( const edm::ParameterSet& );
  /// end-of-job processing
  void endJob();

private:
  /// event processing
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  /// label of source collection
  std::string src_;
  /// verbosity flag
  bool verbose_;
  /// partial statistics
  unsigned long n_, nSum_, n2Sum_;
};

template<typename C>
ObjectCounter<C>::ObjectCounter( const edm::ParameterSet& par ) : 
  src_( par.template getParameter<std::string>( "src" ) ), 
  verbose_( par.template getUntrackedParameter<bool>( "verbose", true ) ),
  n_( 0 ), nSum_( 0 ), n2Sum_( 0 ) {
}

template<typename C>
void ObjectCounter<C>::endJob() {
  double n = 0, n2 = 0, s;
  if ( n_!= 0 ) {
    n = double( nSum_ ) / n_;
    n2 = double ( n2Sum_ ) / n_;
  }
  s = sqrt( n2 - n * n );
  if ( verbose_ ) {
    edm::TypeID id( typeid( typename C::value_type ) );
    std::cout << ">>> collection \"" << src_ << "\" contains (" 
	      << n << " +/- " << s << ") "  
	      << id.friendlyClassName() << " objects" << std::endl;
  }
}

template<typename C>
void ObjectCounter<C>::analyze( const edm::Event& evt, const edm::EventSetup& ) {
  edm::Handle<C> h;
  evt.getByLabel( src_, h );
  if (!h.isValid()) {
    std::cerr << ">>> product: " << src_ << " not found" << std::endl;
  } else {
    int n = h->size();
    nSum_ += n;
    n2Sum_ += ( n * n );
  }
  ++ n_;
}

#endif
