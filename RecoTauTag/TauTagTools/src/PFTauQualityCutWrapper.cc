#include "RecoTauTag/TauTagTools/interface/PFTauQualityCutWrapper.h"

// ****************************************************
// ******   Isolation filters  ************************
// ****************************************************

using namespace reco;

void
PFTauQualityCutWrapper::isolationChargedObjects(const PFTau& pfTau, const Vertex& pv, std::vector<reco::LeafCandidate>& output)
{
   if( isoQCuts.useTracksInsteadOfPF )
   {
      TrackRefVector result = TauTagTools::filteredTracks(pfTau.isolationTracks(), 
                                             isoQCuts.minTrackPt,
                                             isoQCuts.minTrackPixelHits,
                                             isoQCuts.minTrackHits,
                                             isoQCuts.maxTransverseImpactParameter,
                                             isoQCuts.maxTrackChi2,
                                             isoQCuts.maxDeltaZ,
                                             pv,
                                             pv.position().z() ); //????

      size_t nTracks = result.size();
      for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         // this sucks
         int charge         = result[iTrack]->charge();
         math::XYZVector p3 = result[iTrack]->momentum();
         reco::Particle::LorentzVector p4(p3.R(), p3.x(), p3.y(), p3.z());
         output.push_back(reco::LeafCandidate(charge, p4));
      }
   } else
   {
      PFCandidateRefVector result = TauTagTools::filteredPFChargedHadrCands(pfTau.isolationPFChargedHadrCands(), 
                                             isoQCuts.minTrackPt,
                                             isoQCuts.minTrackPixelHits,
                                             isoQCuts.minTrackHits,
                                             isoQCuts.maxTransverseImpactParameter,
                                             isoQCuts.maxTrackChi2,
                                             isoQCuts.maxDeltaZ,
                                             pv,
                                             pv.position().z() ); //????

      size_t nTracks = result.size();
      for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         output.push_back(reco::LeafCandidate(result[iTrack]->charge(), result[iTrack]->p4()));
      }
   }
   
}

void
PFTauQualityCutWrapper::isolationGammaObjects(const PFTau& pfTau, std::vector<reco::LeafCandidate>& output)
{
   PFCandidateRefVector result = TauTagTools::filteredPFGammaCands(pfTau.isolationPFGammaCands(), isoQCuts.minGammaEt);

   size_t nGammas = result.size();
   for(size_t iGamma = 0; iGamma < nGammas; ++iGamma)
   {
      output.push_back(reco::LeafCandidate(result[iGamma]->charge(), result[iGamma]->p4()));
   }
   
}

// ****************************************************
// ******   Signal region filters *********************
// ****************************************************

void
PFTauQualityCutWrapper::signalChargedObjects(const PFTau& pfTau, const Vertex& pv, std::vector<reco::LeafCandidate>& output)
{
   if( signalQCuts.useTracksInsteadOfPF )
   {
      TrackRefVector result = TauTagTools::filteredTracks(pfTau.signalTracks(), 
                                             signalQCuts.minTrackPt,
                                             signalQCuts.minTrackPixelHits,
                                             signalQCuts.minTrackHits,
                                             signalQCuts.maxTransverseImpactParameter,
                                             signalQCuts.maxTrackChi2,
                                             signalQCuts.maxDeltaZ,
                                             pv,
                                             pv.position().z() ); //????

      size_t nTracks = result.size();
      for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         // this sucks
         int charge         = result[iTrack]->charge();
         math::XYZVector p3 = result[iTrack]->momentum();
         reco::Particle::LorentzVector p4(p3.R(), p3.x(), p3.y(), p3.z());
         output.push_back(reco::LeafCandidate(charge, p4));
      }
   } else
   {
      PFCandidateRefVector result = TauTagTools::filteredPFChargedHadrCands(pfTau.signalPFChargedHadrCands(), 
                                             signalQCuts.minTrackPt,
                                             signalQCuts.minTrackPixelHits,
                                             signalQCuts.minTrackHits,
                                             signalQCuts.maxTransverseImpactParameter,
                                             signalQCuts.maxTrackChi2,
                                             signalQCuts.maxDeltaZ,
                                             pv,
                                             pv.position().z() ); //????

      size_t nTracks = result.size();
      for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         output.push_back(reco::LeafCandidate(result[iTrack]->charge(), result[iTrack]->p4()));
      }
   }
   
}

void
PFTauQualityCutWrapper::signalGammaObjects(const PFTau& pfTau, std::vector<reco::LeafCandidate>& output)
{
   PFCandidateRefVector result = TauTagTools::filteredPFGammaCands(pfTau.signalPFGammaCands(), signalQCuts.minGammaEt);

   size_t nGammas = result.size();
   for(size_t iGamma = 0; iGamma < nGammas; ++iGamma)
   {
      output.push_back(reco::LeafCandidate(result[iGamma]->charge(), result[iGamma]->p4()));
   }
   
}
