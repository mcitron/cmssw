//---------------------------------------------
//
//   \class L1MuGMTHWFileReader
//
//   Description: Puts the GMT input information from 
//                a GMT ascii HW testfile into the Event   
//
//
//   $Date: 2006/05/15 13:56:02 $
//   $Revision: 1.1 $
//
//   Author :
//   Tobias Noebauer                 HEPHY Vienna
//   Ivan Mikulec                    HEPHY Vienna
//
//--------------------------------------------------
#ifndef L1TriggerGlobalMuonTrigger_L1MuGMTHWFileReader_h
#define L1TriggerGlobalMuonTrigger_L1MuGMTHWFileReader_h

//---------------
// C++ Headers --
//---------------
#include <fstream>

//----------------------
// Base Class Headers --
//----------------------
#include "FWCore/Sources/interface/ExternalInputSource.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "L1Trigger/GlobalMuonTrigger/src/L1MuGMTInputEvent.h"

//---------------------
//-- Class Interface --
//---------------------

class L1MuGMTHWFileReader : public edm::ExternalInputSource {

 public:
  explicit L1MuGMTHWFileReader(edm::ParameterSet const&,
                               edm::InputSourceDescription const&);

  ~L1MuGMTHWFileReader();

  //read an event from the input stream
  //returns an event with run and event number zero when no more events
  void readNextEvent();

  virtual void setRunAndEventInfo();
  virtual bool produce(edm::Event&);

 private:
   std::ifstream m_in;
   L1MuGMTInputEvent m_evt;
};

#endif // L1TriggerGlobalMuonTrigger_L1MuGMTHWFileReader_h
