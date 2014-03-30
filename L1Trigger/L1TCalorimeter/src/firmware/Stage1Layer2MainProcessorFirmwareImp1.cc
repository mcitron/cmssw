///
/// \class l1t::Stage1Layer2MainProcessorFirmwareImp1
///
///
/// \author: R. Alex Barbieri MIT
///
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2MainProcessorFirmware.h"
#include "CondFormats/L1TObjects/interface/FirmwareVersion.h"

#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2EGammaAlgorithmImp.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2EtSumAlgorithmImp.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2JetAlgorithmImp.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2TauAlgorithmImp.h"

using namespace std;
using namespace l1t;

// Stage1Layer2MainProcessorFirmwareImp1::Stage1Layer2MainProcessorFirmwareImp1(/*const CaloParams & dbPars*/
Stage1Layer2MainProcessorFirmwareImp1::Stage1Layer2MainProcessorFirmwareImp1(const FirmwareVersion & fwv, CaloParams* dbPars) : m_fwv(fwv), m_db(dbPars) {}

Stage1Layer2MainProcessorFirmwareImp1::~Stage1Layer2MainProcessorFirmwareImp1(){};

void Stage1Layer2MainProcessorFirmwareImp1::processEvent(const std::vector<CaloEmCand> & emcands,
						       const std::vector<CaloRegion> & regions,
						       std::vector<EGamma> * egammas,
						       std::vector<Tau> * taus,
						       std::vector<Jet> * jets,
						       std::vector<EtSum> * etsums){

  if (m_fwv.firmwareVersion() == 1)
  { //HI algo
    //m_egAlgo = new Stage1Layer2EGammaAlgorithmImpHI();
    m_sumAlgo = new Stage1Layer2EtSumAlgorithmImpPP(m_db);
    m_jetAlgo = new Stage1Layer2JetAlgorithmImpHI(m_db); //fwv =1 => HI algo
    m_tauAlgo = new Stage1Layer2SingleTrackHI(/*m_db*/);
  }
  else if( m_fwv.firmwareVersion() == 2 )
  { //PP algorithm
    m_egAlgo = new Stage1Layer2EGammaAlgorithmImpPP(/*m_db*/);
    m_sumAlgo = new Stage1Layer2EtSumAlgorithmImpPP(m_db);
    m_jetAlgo = new Stage1Layer2JetAlgorithmImpPP(m_db); //fwv =2 => PP algo
    m_tauAlgo = new Stage1Layer2SingleTrackHI(/*m_db*/); //only for now
  }
  else if( m_fwv.firmwareVersion() == 3 )
  {
    //m_tauAlgo = new Stage1Layer2SingleTrackHI(/*m_db*/);
  }
  else{ // undefined fwv version
    edm::LogError("FWVersionError")
      << "Undefined firmware version passed to Stage1Layer2MainProcessorFirmwareImp1" << std::endl;
    return;
  }

  m_egAlgo->processEvent(emcands, regions, egammas);
  m_sumAlgo->processEvent(regions, etsums);
  m_jetAlgo->processEvent(regions, emcands, jets);
  m_tauAlgo->processEvent(emcands, regions, taus);
}
