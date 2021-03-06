#ifndef SISTRIPCLUSTERIZER_SISTRIPCLUSTERINFO_H
#define SISTRIPCLUSTERIZER_SISTRIPCLUSTERINFO_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include <numeric>

class SiStripNoises;
class SiStripGain;
class SiStripQuality;


class SiStripClusterInfo {

 public:

  SiStripClusterInfo(const SiStripCluster& cluster, 
		     const edm::EventSetup& es, 
		     std::string qualityLabel="");

  const SiStripCluster * cluster() const {return cluster_ptr;}

  uint32_t detId() const      {return cluster()->geographicalId();}
  uint16_t width() const      {return cluster()->amplitudes().size();}
  uint16_t firstStrip() const {return cluster()->firstStrip();}
  float    baryStrip() const  {return cluster()->barycenter();}
  uint16_t maxStrip() const   {return firstStrip() + maxIndex();}
  float    variance() const;

  const std::vector<uint8_t>& stripCharges() const {return cluster()->amplitudes();}
  std::vector<float>          stripGains() const;
  std::vector<float>          stripNoises() const;
  std::vector<float>          stripNoisesRescaledByGain() const;
  std::vector<bool>           stripQualitiesBad() const;

  uint16_t charge() const    {return   accumulate( stripCharges().begin(), stripCharges().end(), uint16_t(0));}
  uint8_t  maxCharge() const {return * max_element(stripCharges().begin(), stripCharges().end());}
  uint16_t maxIndex() const  {return   max_element(stripCharges().begin(), stripCharges().end()) - stripCharges().begin();}
  std::pair<uint16_t,uint16_t> chargeLR() const;
  
  float noise() const               { return calculate_noise(stripNoises());}
  float noiseRescaledByGain() const { return calculate_noise(stripNoisesRescaledByGain());}

  float signalOverNoise() const { return charge()/noiseRescaledByGain(); }

  bool IsAnythingBad() const;
  bool IsApvBad() const;
  bool IsFiberBad() const;
  bool IsModuleBad() const;
  bool IsModuleUsable() const;

  std::vector<SiStripCluster> reclusterize(const edm::ParameterSet&) const;

 private:

  float calculate_noise(const std::vector<float>&) const; 

  const SiStripCluster* cluster_ptr;
  const edm::EventSetup& es; 
  edm::ESHandle<SiStripNoises> noiseHandle;
  edm::ESHandle<SiStripGain> gainHandle;
  edm::ESHandle<SiStripQuality> qualityHandle;
  std::string qualityLabel;

};

#endif
