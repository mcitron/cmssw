#ifndef FastSimulation_ParamL3MuonProducer_ParamL3MuonProducer_h
#define FastSimulation_ParamL3MuonProducer_ParamL3MuonProducer_h

//
// Package:    ParamL3MuonProducer
// Class:      ParamL3MuonProducer
// 
/**\class ParamL3MuonProducer ParamL3MuonProducer.cc FastSimulation/ParamL3MuonProducer/src/ParamL3MuonProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//  Author:  Andrea Perrotta
// Created:  Wed May 02 12:37:24 CET 2007
// $Id: ParamL3MuonProducer.h,v 1.9 2009/03/03 14:06:08 abdullin Exp $
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EDProducer.h"

// FastSimulation headers
class SimpleL1MuGMTCand;
class FML1EfficiencyHandler;
class FML1PtSmearer;
class FML3EfficiencyHandler; 
class FML3PtSmearer;
class FMGLfromL3EfficiencyHandler; 
class FMGLfromL3TKEfficiencyHandler; 
class FMGLfromTKEfficiencyHandler; 

class RandomEngine;

namespace reco { 
  class Muon;
  class Track;
}

namespace edm { 
  class ParameterSet;
  class Event;
  class EventSetup;
}

// Data Formats
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

//
// class declaration
//

class ParamL3MuonProducer : public edm::EDProducer {
   public:

      explicit ParamL3MuonProducer(const edm::ParameterSet&);
      ~ParamL3MuonProducer();

   private:

      int findId(const reco::Track& aTrack) const;

   private:

      const RandomEngine * random;

      typedef std::vector<SimpleL1MuGMTCand*> FML1Muons;
      typedef std::vector<L1MuGMTCand> L1MuonCollection;
      typedef std::vector<l1extra::L1MuonParticle> L1ExtraCollection;

      virtual void beginRun(edm::Run&, edm::EventSetup const& ) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endRun() ;
      void readParameters(const edm::ParameterSet&, const edm::ParameterSet& );
      void reconstruct();
      void loadL1Muons(L1MuonCollection & c, L1ExtraCollection & d) const;
      void loadL3Muons(reco::MuonCollection & c) const;
      void loadGLMuons(reco::MuonCollection & c) const;
    
  // ---------- member data ---------------------------

      FML1Muons  mySimpleL1MuonCands;
      L1ExtraCollection mySimpleL1MuonExtraCands;
      FML1EfficiencyHandler * myL1EfficiencyHandler;
      FML1PtSmearer * myL1PtSmearer;

      reco::MuonCollection  mySimpleL3MuonCands;
      FML1Muons  mySimpleL3MuonSeeds;
      FML3EfficiencyHandler * myL3EfficiencyHandler;
      FML3PtSmearer * myL3PtSmearer;

      reco::MuonCollection  mySimpleGLMuonCands;
      FMGLfromL3EfficiencyHandler * myGLfromL3EfficiencyHandler;
      FMGLfromL3TKEfficiencyHandler * myGLfromL3TKEfficiencyHandler;
      FMGLfromTKEfficiencyHandler * myGLfromTKEfficiencyHandler;
      FML3PtSmearer * myGLPtSmearer;
      
  // ----------- parameters ---------------------------- 
      bool fullPattern_;
      bool doL1_ , doL3_ , doGL_;
      std::string theSimModuleLabel_ , theSimModuleProcess_, theTrkModuleLabel_ ;
      double minEta_ ,  maxEta_;
  // ----------- counters ------------------------------
      int   nMuonTot , nL1MuonTot , nL3MuonTot , nGLMuonTot;

      static double muonMassGeV_ ;

};

#endif
