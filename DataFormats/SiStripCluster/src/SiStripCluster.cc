
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"

SiStripCluster::SiStripCluster( uint32_t detid, const SiStripDigiRange& range) :
  detId_(detid), firstStrip_(range.first->strip())
{

  amplitudes_.reserve( range.second - range.first);
  
  uint16_t lastStrip=0;
  bool firstInloop = true;
  for (SiStripDigiIter i=range.first; i!=range.second; i++) {
    
    /// check if digis consecutive
    if (!firstInloop && i->strip() != lastStrip + 1) {
      for (int j=0; j < i->strip()-(lastStrip+1); j++) {
	amplitudes_.push_back( 0);
      }
    }
    lastStrip = i->strip();
    firstInloop = false;
    
    amplitudes_.push_back(i->adc()); 
  }
}

SiStripCluster::SiStripCluster(const uint32_t& detid, 
			       const uint16_t& firstStrip, 
			       std::vector<uint16_t>::const_iterator begin, 
			       std::vector<uint16_t>::const_iterator end) :

  detId_(detid),
  firstStrip_(firstStrip),
  amplitudes_(begin,end)
{}

float SiStripCluster::barycenter() const{
  int sumx = 0;
  int suma = 0;
  size_t asize = amplitudes_.size();
  for (size_t i=0;i<asize;++i) {
    sumx += (firstStrip_+i)*(amplitudes_[i]);
    suma += amplitudes_[i];
  }
  
  // strip centers are offcet by half pitch w.r.t. strip numbers,
  // so one has to add 0.5 to get the correct barycenter position
  return sumx / static_cast<float>(suma) + 0.5;
}
