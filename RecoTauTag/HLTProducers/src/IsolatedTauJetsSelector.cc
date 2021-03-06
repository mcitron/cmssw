#include "RecoTauTag/HLTProducers/interface/IsolatedTauJetsSelector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Math/GenVector/VectorUtil.h"
//
// class decleration
//


IsolatedTauJetsSelector::IsolatedTauJetsSelector(const edm::ParameterSet& iConfig)
{

  jetSrc = iConfig.getParameter<vtag>("JetSrc");
  pt_min_leadTrack   = iConfig.getParameter<double>("MinimumTransverseMomentumLeadingTrack"); 
useIsolationDiscriminator = iConfig.getParameter<bool>("UseIsolationDiscriminator");
 useInHLTOpen       = iConfig.getParameter<bool>("UseInHLTOpen");
 
  produces<reco::CaloJetCollection>();
  //  produces<reco::IsolatedTauTagInfoCollection>();  
}

IsolatedTauJetsSelector::~IsolatedTauJetsSelector(){ }

void IsolatedTauJetsSelector::produce(edm::Event& iEvent, const edm::EventSetup& iES)
{

  using namespace reco;
  using namespace edm;
  using namespace std;
  
  CaloJetCollection * jetCollectionTmp = new CaloJetCollection;
//    IsolatedTauTagInfoCollection * extendedCollection = new IsolatedTauTagInfoCollection;

  for( vtag::const_iterator s = jetSrc.begin(); s != jetSrc.end(); ++ s ) {
    edm::Handle<IsolatedTauTagInfoCollection> tauJets;
    iEvent.getByLabel( * s, tauJets );
    IsolatedTauTagInfoCollection::const_iterator i = tauJets->begin();
    for(;i !=tauJets->end(); i++ ) {

      if(useInHLTOpen) {
	const CaloJet* pippo = dynamic_cast<const CaloJet*>((i->jet().get()));
	CaloJet* mioPippo = const_cast<CaloJet*>(pippo);
	mioPippo->setPdgId(15);
	if(mioPippo)
	  jetCollectionTmp->push_back(*mioPippo);
	//	extendedCollection->push_back(*(i)); //to  be used later
	//	delete pippo;
      }else{
	const TrackRef leadTk = i->leadingSignalTrack();
	if( !leadTk ) 
	  {}else{
	if(leadTk->pt() >  pt_min_leadTrack) {	   
	  float discriminator = i->discriminator();	  
	  const CaloJet* pippo = dynamic_cast<const CaloJet*>((i->jet().get()));
	  CaloJet* mioPippo = const_cast<CaloJet*>(pippo);
	  mioPippo->setPdgId(15);
	  if(useIsolationDiscriminator && (discriminator > 0) ) {
	    if(mioPippo)
	      jetCollectionTmp->push_back(*mioPippo);
	    // delete pippo;
	  }else if(!useIsolationDiscriminator){
	    if(mioPippo)
	      jetCollectionTmp->push_back(*mioPippo);
	    }
	  }
	}
      }
    }	  
  }


  
  auto_ptr<reco::CaloJetCollection> selectedTaus(jetCollectionTmp);

  iEvent.put(selectedTaus);


}
