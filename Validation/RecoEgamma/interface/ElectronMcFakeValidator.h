
#ifndef Validation_RecoEgamma_ElectronMcFakeValidator_h
#define Validation_RecoEgamma_ElectronMcFakeValidator_h

#include "Validation/RecoEgamma/interface/ElectronValidator.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
class MagneticField;

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

class ElectronMcFakeValidator : public ElectronValidator
 {
  public:

    explicit ElectronMcFakeValidator( const edm::ParameterSet & conf ) ;
    virtual ~ElectronMcFakeValidator() ;

    virtual void beginJob() ;
    virtual void endJob() ;
    virtual void analyze( const edm::Event& e, const edm::EventSetup & c ) ;

  private:

    edm::InputTag electronCollection_;
    edm::InputTag  matchingObjectCollection_;
    edm::InputTag beamSpotTag_;
    bool readAOD_;
    std::string outputFile_ ;

    TrajectoryStateTransform transformer_ ;
    edm::ESHandle<TrackerGeometry> pDD ;
    edm::ESHandle<MagneticField> theMagField ;

    float mcEnergy[10], mcEta[10], mcPhi[10], mcPt[10], mcQ[10] ;
    float superclusterEnergy[10], superclusterEta[10], superclusterPhi[10], superclusterEt[10] ;
    float seedMomentum[10], seedEta[10], seedPhi[10], seedPt[10], seedQ[10] ;

    double maxPt_;
    double maxAbsEta_;
    double deltaR_;

    // histos limits and binning

    int xyz_nbin ;
    int p_nbin ; int p2D_nbin ; double p_max ;
    int pt_nbin ; int pt2D_nbin ; int pteff_nbin ; double pt_max ;
    int fhits_nbin ; double fhits_max ;
    int lhits_nbin ; double lhits_max ;
    int eta_nbin ; int eta2D_nbin ; double eta_min ; double eta_max ;
    int deta_nbin ; double deta_min ; double deta_max ;
    int detamatch_nbin ; int detamatch2D_nbin ; double detamatch_min ; double detamatch_max ;
    int phi_nbin ; int phi2D_nbin ; double phi_min ; double phi_max ;
    int dphi_nbin ; double dphi_min ; double dphi_max ;
    int dphimatch_nbin ; int    dphimatch2D_nbin ; double dphimatch_min ; double dphimatch_max ;
    int eop_nbin ; int eop2D_nbin ; double eop_max ; double eopmaxsht ;
    int mee_nbin ; double mee_min ; double mee_max ;
    int hoe_nbin ; double hoe_min ; double hoe_max ;
    int popmatching_nbin ; double popmatching_min ; double popmatching_max ;

    // histos

    MonitorElement *h1_matchingObjectNum ;
    MonitorElement *h1_recEleNum_;

    MonitorElement *h1_matchingObjectEta;
    MonitorElement *h1_matchingObjectAbsEta;
    MonitorElement *h1_matchingObjectP;
    MonitorElement *h1_matchingObjectPt;
    MonitorElement *h1_matchingObjectPhi;
    MonitorElement *h1_matchingObjectZ;

    MonitorElement *h1_ele_EoverP_all;
    MonitorElement *h1_ele_EseedOP_all;
    MonitorElement *h1_ele_EoPout_all;
    MonitorElement *h1_ele_EeleOPout_all;
    MonitorElement *h1_ele_dEtaSc_propVtx_all;
    MonitorElement *h1_ele_dPhiSc_propVtx_all;
    MonitorElement *h1_ele_dEtaCl_propOut_all;
    MonitorElement *h1_ele_dPhiCl_propOut_all;
    MonitorElement *h1_ele_TIP_all;
    MonitorElement *h1_ele_HoE_all;
    MonitorElement *h1_ele_vertexEta_all;
    MonitorElement *h1_ele_vertexPt_all;
    MonitorElement *h1_ele_mee_all;
    MonitorElement *h1_ele_mee_os;

    MonitorElement *h2_ele_E2mnE1vsMee_all;
    MonitorElement *h2_ele_E2mnE1vsMee_egeg_all;

    MonitorElement *h1_ele_matchingObjectEta_matched;
    MonitorElement *h1_ele_matchingObjectAbsEta_matched;
    MonitorElement *h1_ele_matchingObjectPt_matched;
    MonitorElement *h1_ele_matchingObjectPhi_matched;
    MonitorElement *h1_ele_matchingObjectZ_matched;

    MonitorElement *h1_ele_charge;
    MonitorElement *h2_ele_chargeVsEta;
    MonitorElement *h2_ele_chargeVsPhi;
    MonitorElement *h2_ele_chargeVsPt;
    MonitorElement *h1_ele_vertexP;
    MonitorElement *h1_ele_vertexPt;
    MonitorElement *h2_ele_vertexPtVsEta;
    MonitorElement *h2_ele_vertexPtVsPhi;
    MonitorElement *h1_ele_vertexEta;
    MonitorElement *h2_ele_vertexEtaVsPhi;
    MonitorElement *h1_ele_vertexAbsEta;
    MonitorElement *h1_ele_vertexPhi;
    MonitorElement *h1_ele_vertexX;
    MonitorElement *h1_ele_vertexY;
    MonitorElement *h1_ele_vertexZ;
    MonitorElement *h1_ele_vertexTIP;
    MonitorElement *h2_ele_vertexTIPVsEta;
    MonitorElement *h2_ele_vertexTIPVsPhi;
    MonitorElement *h2_ele_vertexTIPVsPt;

    MonitorElement *h1_ele_PoPmatchingObject;
    MonitorElement *h2_ele_PoPmatchingObjectVsEta;
    MonitorElement *h2_ele_PoPmatchingObjectVsPhi;
    MonitorElement *h2_ele_PoPmatchingObjectVsPt;
    MonitorElement *h1_ele_PoPmatchingObject_barrel;
    MonitorElement *h1_ele_PoPmatchingObject_endcaps;

    MonitorElement *h1_ele_EtaMnEtamatchingObject;
    MonitorElement *h2_ele_EtaMnEtamatchingObjectVsEta;
    MonitorElement *h2_ele_EtaMnEtamatchingObjectVsPhi;
    MonitorElement *h2_ele_EtaMnEtamatchingObjectVsPt;
    MonitorElement *h1_ele_PhiMnPhimatchingObject;
    MonitorElement *h1_ele_PhiMnPhimatchingObject2;
    MonitorElement *h2_ele_PhiMnPhimatchingObjectVsEta;
    MonitorElement *h2_ele_PhiMnPhimatchingObjectVsPhi;
    MonitorElement *h2_ele_PhiMnPhimatchingObjectVsPt;

    MonitorElement *h1_scl_En_ ;
    MonitorElement *h1_scl_EoEmatchingObject_barrel;
    MonitorElement *h1_scl_EoEmatchingObject_endcaps;
    MonitorElement *h1_scl_Et_ ;
    MonitorElement *h2_scl_EtVsEta_ ;
    MonitorElement *h2_scl_EtVsPhi_ ;
    MonitorElement *h2_scl_EtaVsPhi_ ;
    MonitorElement *h1_scl_Eta_ ;
    MonitorElement *h1_scl_Phi_ ;

//    MonitorElement *h1_scl_SigEtaEta_ ;
//    MonitorElement *h1_scl_SigEtaEta_barrel_ ;
//    MonitorElement *h1_scl_SigEtaEta_endcaps_ ;
    MonitorElement *h1_scl_SigIEtaIEta_ ;
    MonitorElement *h1_scl_SigIEtaIEta_barrel_ ;
    MonitorElement *h1_scl_SigIEtaIEta_endcaps_ ;
    MonitorElement *h1_scl_E1x5_ ;
    MonitorElement *h1_scl_E1x5_barrel_ ;
    MonitorElement *h1_scl_E1x5_endcaps_ ;
    MonitorElement *h1_scl_E2x5max_ ;
    MonitorElement *h1_scl_E2x5max_barrel_ ;
    MonitorElement *h1_scl_E2x5max_endcaps_ ;
    MonitorElement *h1_scl_E5x5_ ;
    MonitorElement *h1_scl_E5x5_barrel_ ;
    MonitorElement *h1_scl_E5x5_endcaps_ ;

    MonitorElement *h1_ele_ambiguousTracks;
    MonitorElement *h2_ele_ambiguousTracksVsEta;
    MonitorElement *h2_ele_ambiguousTracksVsPhi;
    MonitorElement *h2_ele_ambiguousTracksVsPt;
    MonitorElement *h1_ele_foundHits;
    MonitorElement *h1_ele_foundHits_barrel;
    MonitorElement *h1_ele_foundHits_endcaps;
    MonitorElement *h2_ele_foundHitsVsEta;
    MonitorElement *h2_ele_foundHitsVsPhi;
    MonitorElement *h2_ele_foundHitsVsPt;
    MonitorElement *h1_ele_lostHits;
    MonitorElement *h1_ele_lostHits_barrel;
    MonitorElement *h1_ele_lostHits_endcaps;
    MonitorElement *h2_ele_lostHitsVsEta;
    MonitorElement *h2_ele_lostHitsVsPhi;
    MonitorElement *h2_ele_lostHitsVsPt;
    MonitorElement *h1_ele_chi2;
    MonitorElement *h1_ele_chi2_barrel;
    MonitorElement *h1_ele_chi2_endcaps;
    MonitorElement *h2_ele_chi2VsEta;
    MonitorElement *h2_ele_chi2VsPhi;
    MonitorElement *h2_ele_chi2VsPt;

    MonitorElement *h1_ele_PinMnPout;
    MonitorElement *h1_ele_PinMnPout_mode;
    MonitorElement *h2_ele_PinMnPoutVsEta_mode;
    MonitorElement *h2_ele_PinMnPoutVsPhi_mode;
    MonitorElement *h2_ele_PinMnPoutVsPt_mode;
    MonitorElement *h2_ele_PinMnPoutVsE_mode;
    MonitorElement *h2_ele_PinMnPoutVsChi2_mode;

    MonitorElement *h1_ele_outerP;
    MonitorElement *h1_ele_outerP_mode;
    MonitorElement *h2_ele_outerPVsEta_mode;
    MonitorElement *h1_ele_outerPt;
    MonitorElement *h1_ele_outerPt_mode;
    MonitorElement *h2_ele_outerPtVsEta_mode;
    MonitorElement *h2_ele_outerPtVsPhi_mode;
    MonitorElement *h2_ele_outerPtVsPt_mode;
    MonitorElement *h1_ele_EoP;
    MonitorElement *h1_ele_EoP_barrel;
    MonitorElement *h1_ele_EoP_endcaps;
    MonitorElement *h1_ele_EoP_eg;
    MonitorElement *h1_ele_EoP_eg_barrel;
    MonitorElement *h1_ele_EoP_eg_endcaps;
    MonitorElement *h2_ele_EoPVsEta;
    MonitorElement *h2_ele_EoPVsPhi;
    MonitorElement *h2_ele_EoPVsE;
    MonitorElement *h1_ele_EseedOP;
    MonitorElement *h1_ele_EseedOP_barrel;
    MonitorElement *h1_ele_EseedOP_endcaps;
    MonitorElement *h1_ele_EseedOP_eg;
    MonitorElement *h1_ele_EseedOP_eg_barrel;
    MonitorElement *h1_ele_EseedOP_eg_endcaps;
    MonitorElement *h2_ele_EseedOPVsEta;
    MonitorElement *h2_ele_EseedOPVsPhi;
    MonitorElement *h2_ele_EseedOPVsE;
    MonitorElement *h1_ele_EoPout;
    MonitorElement *h1_ele_EoPout_barrel;
    MonitorElement *h1_ele_EoPout_endcaps;
    MonitorElement *h1_ele_EoPout_eg;
    MonitorElement *h1_ele_EoPout_eg_barrel;
    MonitorElement *h1_ele_EoPout_eg_endcaps;
    MonitorElement *h2_ele_EoPoutVsEta;
    MonitorElement *h2_ele_EoPoutVsPhi;
    MonitorElement *h2_ele_EoPoutVsE;
    MonitorElement *h1_ele_EeleOPout;
    MonitorElement *h1_ele_EeleOPout_barrel;
    MonitorElement *h1_ele_EeleOPout_endcaps;
    MonitorElement *h1_ele_EeleOPout_eg;
    MonitorElement *h1_ele_EeleOPout_eg_barrel;
    MonitorElement *h1_ele_EeleOPout_eg_endcaps;
    MonitorElement *h2_ele_EeleOPoutVsEta;
    MonitorElement *h2_ele_EeleOPoutVsPhi;
    MonitorElement *h2_ele_EeleOPoutVsE;

    MonitorElement *h1_ele_dEtaSc_propVtx;
    MonitorElement *h1_ele_dEtaSc_propVtx_barrel;
    MonitorElement *h1_ele_dEtaSc_propVtx_endcaps;
    MonitorElement *h1_ele_dEtaSc_propVtx_eg;
    MonitorElement *h1_ele_dEtaSc_propVtx_eg_barrel;
    MonitorElement *h1_ele_dEtaSc_propVtx_eg_endcaps;
    MonitorElement *h2_ele_dEtaScVsEta_propVtx;
    MonitorElement *h2_ele_dEtaScVsPhi_propVtx;
    MonitorElement *h2_ele_dEtaScVsPt_propVtx;
    MonitorElement *h1_ele_dPhiSc_propVtx;
    MonitorElement *h1_ele_dPhiSc_propVtx_barrel;
    MonitorElement *h1_ele_dPhiSc_propVtx_endcaps;
    MonitorElement *h1_ele_dPhiSc_propVtx_eg;
    MonitorElement *h1_ele_dPhiSc_propVtx_eg_barrel;
    MonitorElement *h1_ele_dPhiSc_propVtx_eg_endcaps;
    MonitorElement *h2_ele_dPhiScVsEta_propVtx;
    MonitorElement *h2_ele_dPhiScVsPhi_propVtx;
    MonitorElement *h2_ele_dPhiScVsPt_propVtx;
    MonitorElement *h1_ele_dEtaCl_propOut;
    MonitorElement *h1_ele_dEtaCl_propOut_barrel;
    MonitorElement *h1_ele_dEtaCl_propOut_endcaps;
    MonitorElement *h1_ele_dEtaCl_propOut_eg;
    MonitorElement *h1_ele_dEtaCl_propOut_eg_barrel;
    MonitorElement *h1_ele_dEtaCl_propOut_eg_endcaps;
    MonitorElement *h2_ele_dEtaClVsEta_propOut;
    MonitorElement *h2_ele_dEtaClVsPhi_propOut;
    MonitorElement *h2_ele_dEtaClVsPt_propOut;
    MonitorElement *h1_ele_dPhiCl_propOut;
    MonitorElement *h1_ele_dPhiCl_propOut_barrel;
    MonitorElement *h1_ele_dPhiCl_propOut_endcaps;
    MonitorElement *h1_ele_dPhiCl_propOut_eg;
    MonitorElement *h1_ele_dPhiCl_propOut_eg_barrel;
    MonitorElement *h1_ele_dPhiCl_propOut_eg_endcaps;
    MonitorElement *h2_ele_dPhiClVsEta_propOut;
    MonitorElement *h2_ele_dPhiClVsPhi_propOut;
    MonitorElement *h2_ele_dPhiClVsPt_propOut;
    MonitorElement *h1_ele_dEtaEleCl_propOut;
    MonitorElement *h1_ele_dEtaEleCl_propOut_barrel;
    MonitorElement *h1_ele_dEtaEleCl_propOut_endcaps;
    MonitorElement *h1_ele_dEtaEleCl_propOut_eg;
    MonitorElement *h1_ele_dEtaEleCl_propOut_eg_barrel;
    MonitorElement *h1_ele_dEtaEleCl_propOut_eg_endcaps;
    MonitorElement *h2_ele_dEtaEleClVsEta_propOut;
    MonitorElement *h2_ele_dEtaEleClVsPhi_propOut;
    MonitorElement *h2_ele_dEtaEleClVsPt_propOut;
    MonitorElement *h1_ele_dPhiEleCl_propOut;
    MonitorElement *h1_ele_dPhiEleCl_propOut_barrel;
    MonitorElement *h1_ele_dPhiEleCl_propOut_endcaps;
    MonitorElement *h1_ele_dPhiEleCl_propOut_eg;
    MonitorElement *h1_ele_dPhiEleCl_propOut_eg_barrel;
    MonitorElement *h1_ele_dPhiEleCl_propOut_eg_endcaps;
    MonitorElement *h2_ele_dPhiEleClVsEta_propOut;
    MonitorElement *h2_ele_dPhiEleClVsPhi_propOut;
    MonitorElement *h2_ele_dPhiEleClVsPt_propOut;

    MonitorElement *h1_ele_seed_dphi2_;
    MonitorElement *h2_ele_seed_dphi2VsEta_;
    MonitorElement *h2_ele_seed_dphi2VsPt_ ;
    MonitorElement *h1_ele_seed_drz2_;
    MonitorElement *h2_ele_seed_drz2VsEta_;
    MonitorElement *h2_ele_seed_drz2VsPt_;
    MonitorElement *h1_ele_seed_subdet2_;

    MonitorElement *h1_ele_classes;
    MonitorElement *h1_ele_eta;
    MonitorElement *h1_ele_eta_golden;
    MonitorElement *h1_ele_eta_bbrem;
    MonitorElement *h1_ele_eta_narrow;
    MonitorElement *h1_ele_eta_shower;

    MonitorElement *h1_ele_HoE;
    MonitorElement *h1_ele_HoE_barrel;
    MonitorElement *h1_ele_HoE_endcaps;
    MonitorElement *h1_ele_HoE_eg;
    MonitorElement *h1_ele_HoE_eg_barrel;
    MonitorElement *h1_ele_HoE_eg_endcaps;
    MonitorElement *h1_ele_HoE_fiducial;
    MonitorElement *h2_ele_HoEVsEta;
    MonitorElement *h2_ele_HoEVsPhi;
    MonitorElement *h2_ele_HoEVsE;

    MonitorElement *h1_ele_fbrem;
    MonitorElement *p1_ele_fbremVsEta_mode;
    MonitorElement *p1_ele_fbremVsEta_mean;

    MonitorElement *h2_ele_PinVsPoutGolden_mode;
    MonitorElement *h2_ele_PinVsPoutShowering_mode;
    MonitorElement *h2_ele_PinVsPoutGolden_mean;
    MonitorElement *h2_ele_PinVsPoutShowering_mean;
    MonitorElement *h2_ele_PtinVsPtoutGolden_mode;
    MonitorElement *h2_ele_PtinVsPtoutShowering_mode;
    MonitorElement *h2_ele_PtinVsPtoutGolden_mean;
    MonitorElement *h2_ele_PtinVsPtoutShowering_mean;
    MonitorElement *h1_scl_EoEmatchingObjectGolden_barrel;
    MonitorElement *h1_scl_EoEmatchingObjectGolden_endcaps;
    MonitorElement *h1_scl_EoEmatchingObjectShowering_barrel;
    MonitorElement *h1_scl_EoEmatchingObjectShowering_endcaps;

    MonitorElement *h1_ele_mva;
    MonitorElement *h1_ele_provenance;

    MonitorElement *h1_ele_tkSumPt_dr03;
    MonitorElement *h1_ele_ecalRecHitSumEt_dr03;
    MonitorElement *h1_ele_hcalTowerSumEt_dr03_depth1;
    MonitorElement *h1_ele_hcalTowerSumEt_dr03_depth2;
    MonitorElement *h1_ele_tkSumPt_dr04;
    MonitorElement *h1_ele_ecalRecHitSumEt_dr04;
    MonitorElement *h1_ele_hcalTowerSumEt_dr04_depth1;
    MonitorElement *h1_ele_hcalTowerSumEt_dr04_depth2;

 } ;

#endif



